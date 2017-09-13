import sys
import os

from .pyaccessories.SaveLoad import SaveLoad
from .Utilities import Values, Keys, FileExtension
from .Encryption import Encryption


class Setup(object):

    def __init__(self, time_log, custom_terms):
        """
        :param time_log: log the can be written to
        :param custom_terms: dictionary of custom term information with format - 
        Key: used to index the value to the config file for setup
        Value: 3 Item Tuple ("default value", ask user" - i.e. True/False, "type of value" - i.e. str, int....)
        - A value of None is the default for all parts except for "Ask" which is True
        """
        # save the time log
        self.timelog = time_log

        script_dir = sys.path[0]
        # create the path to the json config file
        self.config_json = os.path.join(script_dir, FileExtension.config_json)
        # create a save load object to write and read to the created config file
        self.config_loader = SaveLoad(self.config_json, create=True)

        # Create and setup a tuple dictionary that has the default configuration terms for processing
        self.default_terms = dict()
        self.setup_default_terms()

        # Define and set the final values from the user for all default fields
        self.default_config_values = self.get_config_values(self.default_terms)

        # Define and set the final values from the user for all custom fields
        self.custom_config_values = self.get_config_values(custom_terms)

        # Get the api_key and first_run once the values for the program have been inputted by the user
        self.api_key = self.default_config_values[Keys.redmine_api_key]
        self.nas_mnt = self.default_config_values[Keys.nas_mount]
        self.seconds_between_check = self.default_config_values[Keys.secs_between_checks]
        self.first_run = self.default_config_values[Keys.first_run]

    def setup_default_terms(self):
        """
        Set the default terms to run the Redmine API to a terms dictionary
        """
        # Setting config with pattern -
        # default_dict["json config key"] = ("Default Value", "Ask User", "Value Type")

        self.default_terms[Keys.first_run] = (Values.first_run, False, None)
        self.default_terms[Keys.nas_mount] = (Values.nas_mount_path, True, str)
        self.default_terms[Keys.secs_between_checks] = (Values.check_time, True, int)
        self.default_terms[Keys.redmine_api_key] = ('none', False, str)

    def get_config_values(self, default_terms):
        """Read current values from the configuration file, if there is no configuration file create one and write
        either the default calue to the config file or ask for user input"""
        config_values = dict()

        for key in default_terms:
            config_values[key] = self.config_loader.get(key, default=default_terms[key][0], ask=default_terms[key][1],
                                                        get_type=default_terms[key][2])
        return config_values

    def set_api_key(self, force):
        """
        Sets the API key for the automation process either from the config file or the users input
        :param force: determines whether or not the user is asked for the RedmineAPI key or not
        """

        if self.first_run == 'yes' and force:
            raise ValueError('Need redmine API key!')
        elif self.first_run == 'yes':
            input_api_key = 'y'
        elif not self.first_run == 'yes' and force:
            input_api_key = 'n'
        else:
            self.timelog.time_print("Would you like to set the redmine api key? (y/n)")
            input_api_key = input()

        if input_api_key == 'y':
            self.timelog.time_print("Enter your redmine api key (will be encrypted to file)")
            self.api_key = input()
            # Encode and send to json file
            self.config_loader.redmine_api_key_encrypted = Encryption.encode(Values.encryption_key,
                                                                             self.api_key).decode('utf-8')
            self.config_loader.first_run = 'no'
            self.config_loader.dump(self.config_json)
        else:
            # Import and decode from file
            self.timelog.time_print("Used Redmine API key from the json file.")
            self.api_key = Encryption.decode(Values.encryption_key, self.api_key)

        import re
        if not re.match(r'^[a-z0-9]{40}$', self.api_key):
            self.timelog.time_print("Invalid Redmine API key!")
            exit(1)

    def get_custom_term_values(self):
        """Returns a dictionary of the custom terms put in the config file"""
        return self.custom_config_values
