# SNVPhyl Automator

This set of scripts allow for automation of SNVPhyl runs via Redmine by submitting the jobs
through SLURM.

To run:
- Go to the head node's home folder
- `git clone --recursive this_repository`
- add your Redmine API key to upload_file.py where it says INSERT API KEY HERE.
- `source /mnt/nas/Virtual_Environments/Generic_Redmine/bin/activate`
- `python SNVPhylSlurm_Run.py`
- You'll be asked for a bunch of parameters. You should be safe hitting enter and leaving them at defaults.
- When asked for your api key, enter it. (Found on redmine under 'My Account')