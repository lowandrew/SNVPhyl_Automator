from . pyaccessories.SaveLoad import SaveLoad
from . RedmineAPI import RedmineInterface
from . Utilities import Values, FileExtension
from . RedmineIssue import Issue
import os


class RedmineAccess:

    def __init__(self, timelog, api_key):
        """
        Sets the redmine api key to access redmine
        :param api_key: 
        """
        self.timelog = timelog
        self.api_key = api_key
        self.redmine_api = RedmineInterface('http://redmine.biodiversity.agr.gc.ca/', self.api_key)

        import sys
        script_dir = sys.path[0]  # copy current path

        # creates a loader object to write to the json file of responded issues
        self.issue_loader = SaveLoad(os.path.join(script_dir, FileExtension.issues_json),
                                     create=True)
        # creates a list of already responded to issue to query before starting new tasks
        self.rm_responded_issues = set(self.issue_loader.get(Values.responded_issues, default=[], ask=False))

    def retrieve_issues(self, issue_status, issue_title):
        """
        Retrieve issues from Redmine specified with a specific issue status and title 
        :param issue_status: Status of the issue - must match the text as seen on Redmine
        :param issue_title: Title of the issue - must be entered in all lower case
        :return: The issues associated with the above
        """
        self.timelog.time_print("Checking for requests matching Issue Status: %s and Issue Title: %s" % (issue_status,
                                                                                                         issue_title))
        data = self.redmine_api.get_new_issues('cfia')
        found_issues = []

        # find all 'issues' on redmine, add them to data
        # Sort through all the issues with status that has been specified and add them to found
        for issue in data['issues']:
            if issue['id'] not in self.rm_responded_issues and issue['status']['name'] == issue_status:
                if issue['subject'].lower().rstrip() == issue_title:
                    found_issues.append(Issue(issue))

        self.timelog.time_print("Found %d new issue(s)..." % len(found_issues))

        # returns number of issues
        return found_issues

    def get_specified_attachment_types(self, issue, extn='.txt', decode=True):
        """
        :param issue: the issue the files are attached on
        :param extn: the extension type of the files you would like returned
        :param decode: whether or not to decode the file once downloaded
        :return: a list of all the attachments of a Redmine issue with the specified extension
        """
        # create a list of all attachments
        attachments = self.redmine_api.get_issue_data(issue.id)['issue']['attachments']
        files = list()
        for attachment in attachments:
            if attachment['filename'].endswith(extn):
                self.timelog.time_print('Found the attachment %s, downloading to the list...' % attachment['filename'])
                files.append(self.redmine_api.download_file(attachment['content_url'], decode))

        return files

    def get_attached_text_file(self, issue, index):
        """
        Return an attached text file to the user - includes .txt , .csv and .tsv files      
        :param issue: the issue the files are attached on
        :param index: the specific file on the Redmine Issue you would like to receive i.e. [0] = 1 
        :return: return the downloaded file or none if no file is attached
        """

        try:
            # create a list of all attachments
            attachments = self.redmine_api.get_issue_data(issue.id)['issue']['attachments']
            if len(attachments) < 0:
                return None
        except KeyError:
            # If error due to no attachments, return none
            return None

        # could be another layer to the extension if needed or in future for specific types
        # if attachment['filename'].endswith('.txt'):

        try:
            # try to get the file name of the corresponding index, if it does not exist return nothing
            file_name = attachments[index]['filename']
        except IndexError:
            return None

        # Log the file being downloaded
        self.timelog.time_print("Found the attachment to the Redmine Request: %s" % file_name)
        self.timelog.time_print("Downloading file.....")

        file = self.redmine_api.download_file(attachments[index]['content_url'])
        return file

    def get_attached_files(self, issue):
        """
        Used to return all attached files and run custom commands on the returned list
        :param issue: takes an issue object passed through
        :return: returns a list/dict of all the attachments associated with the selected issue
        """
        redmine_data = self.redmine_api.get_issue_data(issue.id)
        return redmine_data['issue']['attachments']

    def log_new_issue(self, issue):
        """
        Add files to the responded log so the action will not be performed again
        """
        self.rm_responded_issues.add(issue.id)
        self.issue_loader.responded_issues = list(self.rm_responded_issues)
        self.issue_loader.dump()

    def update_status_inprogress(self, issue, message=""):
        """
        Updates the issue to In Progress, while posting a specified message
        """
        self.redmine_api.update_issue(issue.id, notes=issue.redmine_msg + message, status_change=2)

    def update_issue_to_author(self, issue, message=""):
        """
        Updates the issue to Feedback and assigns the Redmine Request back to the Author while updating with a message
        """
        self.redmine_api.update_issue(issue.id, notes=issue.redmine_msg + message, status_change=4,
                                      assign_to_id=issue.author_id)
