from RedmineAPI import Access
from RedmineAPI import Utilities
from RedmineAPI import RedmineAPI

import sys
import shutil

# Setup redmine access so files can be uploaded.
timelog = Utilities.create_time_log()  # I think this is necessary
redmine = Access.RedmineAccess(timelog, 'INSERT API KEY HERE')  # Un-hardcode this at some point.

# Zip the output file.
shutil.make_archive('SNVPhyl_' + sys.argv[1], 'zip', 'output')
# Upload the file to Redmine.
redmine.redmine_api.upload_file('SNVPhyl_' + sys.argv[1] + '.zip', int(sys.argv[1]), 'application/zip',
                                file_name_once_uploaded='SNVPhyl_' + sys.argv[1] + '.zip')