from RedmineAPI.Utilities import FileExtension, create_time_log
import shutil
import os
from RedmineAPI.Access import RedmineAccess
from RedmineAPI.Configuration import Setup
import glob
from biotools import mash

from Utilities import CustomKeys, CustomValues


def check_distances(ref_fasta, fastq_folder):
    bad_fastqs = list()
    fastqs = glob.glob(os.path.join(fastq_folder, '*R1*'))
    for fastq in fastqs:
        mash.dist(fastq, ref_fasta)
        mash_output = mash.read_mash_output('distances.tab')
        print(mash_output[0].reference, mash_output[0].query, str(mash_output[0].distance))
        if mash_output[0].distance > 0.06:  # May need to adjust this value.
            bad_fastqs.append(mash_output[0].reference)
    return bad_fastqs


class Automate(object):

    def __init__(self, force):

        # create a log, can be written to as the process continues
        self.timelog = create_time_log(FileExtension.runner_log)

        # Key: used to index the value to the config file for setup
        # Value: 3 Item Tuple ("default value", ask user" - i.e. True/False, "type of value" - i.e. str, int....)
        # A value of None is the default for all parts except for "Ask" which is True
        #custom_terms = {CustomKeys.key_name: (CustomValues.value_name, True, str)}  # *** can be more than 1 ***
        custom_terms = dict()
        # Create a RedmineAPI setup object to create/read/write to the config file and get default arguments
        setup = Setup(time_log=self.timelog, custom_terms=custom_terms)
        setup.set_api_key(force)

        # Custom terms saved to the config after getting user input
        # self.custom_values = setup.get_custom_term_values()
        # *** can be multiple custom values variable, just use the key from above to reference the inputted value ***
        # self.your_custom_value_name = self.custom_values[CustomKeys.key_name]

        # Default terms saved to the config after getting user input
        self.seconds_between_checks = setup.seconds_between_check
        self.nas_mnt = setup.nas_mnt
        self.redmine_api_key = setup.api_key

        # Initialize Redmine wrapper
        self.access_redmine = RedmineAccess(self.timelog, self.redmine_api_key)

        self.botmsg = '\n\n_I am a bot. This action was performed automatically._'  # sets bot message
        # Subject name and Status to be searched on Redmine
        self.issue_title = 'snvphyl'  # must be a lower case string to validate properly
        self.issue_status = 'New'

    def timed_retrieve(self):
        """
        Continuously search Redmine in intervals for the inputted period of time, 
        Log errors to the log file as they occur
        """
        import time
        while True:
            # Get issues matching the issue status and subject
            found_issues = self.access_redmine.retrieve_issues(self.issue_status, self.issue_title)
            # Respond to the issues in the list 1 at a time
            while len(found_issues) > 0:
                self.respond_to_issue(found_issues.pop(len(found_issues) - 1))
            self.timelog.time_print("Waiting for the next check.")
            time.sleep(self.seconds_between_checks)

    def respond_to_issue(self, issue):
        """
        Run the desired automation process on the inputted issue, if there is an error update the author
        :param issue: Specified Redmine issue information
        """
        self.timelog.time_print("Found a request to run. Subject: %s. ID: %s" % (issue.subject, str(issue.id)))
        self.timelog.time_print("Adding to the list of responded to requests.")
        self.access_redmine.log_new_issue(issue)

        try:
            issue.redmine_msg = "Beginning the process for: %s" % issue.subject
            self.access_redmine.update_status_inprogress(issue, self.botmsg)

            ##########################################################################################
            # Make the bio_request folder
            os.makedirs('/mnt/nas/bio_requests/' + str(issue.id))
            # Remember the directory we're in.
            work_dir = '/mnt/nas/bio_requests/' + str(issue.id)
            current_dir = os.getcwd()
            des = issue.description.split('\n')
            # Make our fastq directory.
            os.makedirs(work_dir + '/fastqs')
            compare = False
            # Iterate through description to try to figure out reference/compare seqIDs.
            queries = list()
            reference = list()
            for item in des:
                item = item.upper()
                if 'COMPARE' in item:
                    compare = True
                if compare:
                    queries.append(item.replace('\r', ''))
                else:
                    reference.append(item.replace('\r', ''))
            # Extract reference fasta. Should only be one, but if someone messes up and does more than one, the
            # first one that glob finds will be treated as the reference.
            f = open(work_dir + '/seqid.txt', 'w')
            for item in reference:
                f.write(item + '\n')
            f.close()
            cmd = 'python2 /mnt/nas/WGSspades/file_extractor.py {}/seqid.txt {} /mnt/nas/'.format(work_dir, work_dir)
            os.system(cmd)
            # Extract query fastqs. Need to do in both MiSeq Backup and External MiSeq Backup, as those scripts don't
            # look everywhere.
            f = open(work_dir + '/seqid.txt', 'w')
            for item in queries:
                f.write(item + '\n')
            f.close()
            os.chdir('/mnt/nas/MiSeq_Backup')
            cmd = 'python2 /mnt/nas/MiSeq_Backup/file_extractor.py {}/seqid.txt {} '.format(work_dir, work_dir + '/fastqs')
            os.system(cmd)
            os.chdir('/mnt/nas/External_MiSeq_Backup')
            cmd = 'python2 /mnt/nas/External_MiSeq_Backup/file_extractor.py {}/seqid.txt {} '.format(work_dir, work_dir + '/fastqs')
            os.system(cmd)
            # Before we get going, do some MASHing to make sure that all the files are close to the reference.
            # In the event that some files aren't, list them?
            ref = glob.glob(work_dir + '/*.fasta')[0]
            bad_fastqs = check_distances(ref, os.path.join(work_dir, 'fastqs'))
            if bad_fastqs:
                outstr = ''
                for fastq in bad_fastqs:
                    fastq = os.path.split(fastq)[-1].split('_')[0]
                    outstr += fastq + '\n'
                self.access_redmine.update_issue_to_author(issue, '\nWARNING: MASH screening suggests that the following'
                                                                  ' samples may be too far from the reference. '
                                                                  'You may want to create a new issue without them'
                                                                  ' and try again. Samples: {}'.format(outstr))
            # Get back to where we were.
            os.chdir(current_dir)
            # Now set up the snvphyl batch script and submit it.
            f = open('snvphyl.sh')
            lines = f.readlines()
            f.close()
            # Modify the batch script. Print everything to where it needs to be.
            f = open(work_dir + '/' + str(issue.id) + '.sh', 'w')
            # Figure out how many processors to use: one per sample, max of 48, because that's our smallest compute node.
            if len(queries) > 48:
                cpu_request = 48
            else:
                cpu_request = len(queries)
            # Cap the memory request at 192 gigs, since that's the memory we have available.
            if len(queries) > 19:
                mem_request = 19
            else:
                mem_request = len(queries)
            for line in lines:
                if 'job_%j' in line:
                    line = line.replace('job', 'biorequest_' + str(issue.id) + '_job')
                elif 'ntasks' in line:
                    line = line.replace('14', str(cpu_request))
                elif 'mem' in line:
                    line = line.replace('40000', str(mem_request*10000))
                f.write(line)
            # Glob the bio_request folder for fastas to use as reference, and take the first one.
            ref = glob.glob(work_dir + '/*.fasta')[0]
            # SNVPhyl command.
            f.write('python /mnt/nas/slurmtest/snvphyl-galaxy-cli/bin/snvphyl.py --deploy-docker --fastq-dir {} '
                    '--reference-file {} --min-coverage 5 --output-dir {} --docker-port {} --docker-cpus {}\n'.format(work_dir + '/fastqs', ref,
                                                                                   work_dir + '/output', str(issue.id), str(cpu_request)))
            f.write('cd /mnt/nas/bio_requests/' + str(issue.id) + '\n')
            f.write('deactivate\n')
            f.write('source /mnt/nas/Virtual_Environments/Generic_Redmine/bin/activate\n')
            f.write('python3 upload_file.py {} \n'.format(str(issue.id)))
            f.write('rm -rf upload_file.py fastqs *fasta RedmineAPI running_logs *json *.txt\n')
            f.close()
            # Copy things to the bio_request folder that are needed to upload the result file to redmine.
            # They'll be cleaned up by the batch script.
            shutil.copy('upload_file.py', work_dir + '/upload_file.py')
            shutil.copytree('RedmineAPI', work_dir + '/RedmineAPI')
            # Submit the batch script to slurm.
            cmd = 'sbatch {}'.format(work_dir + '/' + str(issue.id) + '.sh')
            os.system(cmd)
            ##########################################################################################
            self.completed_response(issue)

        except Exception as e:
            import traceback
            self.timelog.time_print("[Warning] The automation process had a problem, continuing redmine api anyways.")
            self.timelog.time_print("[Automation Error Dump]\n" + traceback.format_exc())
            # Send response
            issue.redmine_msg = "There was a problem with your request. Please create a new issue on" \
                                " Redmine to re-run it.\n%s" % traceback.format_exc()
            # Set it to feedback and assign it back to the author
            self.access_redmine.update_issue_to_author(issue, self.botmsg)

    def completed_response(self, issue):
        """
        Update the issue back to the author once the process has finished
        :param issue: Specified Redmine issue the process has been completed on
        """
        # Assign the issue back to the Author
        self.timelog.time_print("Assigning the issue: %s back to the author." % str(issue.id))

        issue.redmine_msg = "Your job has been submitted to the OLC Compute cluster. This issue will be updated" \
                            " with results once the job is complete."
        # Update author on Redmine
        self.access_redmine.update_issue_to_author(issue, self.botmsg)

        # Log the completion of the issue including the message sent to the author
        self.timelog.time_print("\nMessage to author - %s\n" % issue.redmine_msg)
        self.timelog.time_print("Completed Response to issue %s." % str(issue.id))
        self.timelog.time_print("The next request will be processed once available")
