# Written by Devon Mack April 2017
# This class can be used to easily store variables to a file and load them back later
# then it asks the user to input it and updates the file
import json


class SaveLoad(object):
    def __init__(self, file_name=None, create=False):
        if file_name is not None:
            self.load(file_name, create)

        # Save the file name for later use
        self.file_name_saved = file_name

    def dump(self, file_name=None):
        """Dumps this object into file as JSON."""
        # Dump the json
        f = open(self.__get_saved_filename(file_name), "w")
        file_name_saved = self.__dict__.pop('file_name_saved', None)  # So that it doesn't get saved to file
        json.dump(self.__dict__, f, sort_keys=True, indent=4, separators=(',', ': '))
        self.file_name_saved = file_name_saved
        f.close()

    def load(self, file_name=None, create=False):
        """Loads JSON from file into this object."""
        u_file_name = self.__get_saved_filename(file_name)
        try:
            f = open(u_file_name, "r")
            try:
                self.__dict__ = json.load(f)
            except json.decoder.JSONDecodeError as e:
                # If it's an empty file
                if "Expecting value: " in repr(e) and "line 1 column 1 (char 0)" in repr(e):
                    pass
                else:
                    print("Invalid JSON!")
                    raise
            f.close()
            # Successfully loaded file
            return True
        except FileNotFoundError:
            if create:
                print("[SaveLoad] Creating file " + u_file_name)
                f = open(u_file_name, "w")
                f.close()
                # File doesn't exist
                return False
            else:
                raise

    def get(self, variable, file_name=None, default=None, ask=True, get_type=None):
        """Returns the value of the variable specified. If the variable doesn't exist then it will ask the user to input
        this variable. It will then dump the variable to file_name, or the last file_name used."""
        if get_type is not int and get_type is not float and get_type is not str and get_type is not None:
            raise ValueError("get_type must be int, float, or str.")

        if variable in self.__dict__:
            pass
        elif ask:
            u_file_name = self.__get_saved_filename(file_name)
            if default is None:
                self.__dict__[variable] = input("%s not in %s, please enter it here:\n" % (variable, u_file_name))
            else:
                self.__dict__[variable] = input("%s not in %s, please enter it here (press enter for default option"
                                                " \"%s\"):\n" % (variable, u_file_name, default))
                if self.__dict__[variable] == "":
                    self.__dict__[variable] = default
            if get_type is not None:
                before = self.__dict__[variable]
                if get_type is int:
                    after = int(before)
                elif get_type is float:
                    after = float(before)
                elif get_type is str:
                    after = str(before)
                else:
                    raise ValueError('This shouldn\'t be running')
                self.__dict__[variable] = after

            self.dump()
            print("Updated %s" % u_file_name)
        else:
            u_file_name = self.__get_saved_filename(file_name)
            if default is not None:
                self.__dict__[variable] = default
            else:
                raise ValueError("Missing variable %s in %s!" % (variable, u_file_name))

            self.dump()
            print("Updated %s" % u_file_name)

        return self.__dict__[variable]

    def __get_saved_filename(self, file_name):
        """If the file name inputted is None then this will return whatever is saved. If nothing is saved then it will
        throw an error. If the file name inputted is not none it will return what is inputted back."""
        if file_name is not None:  # They entered a file name so save it
            return file_name
        else:  # No file name entered so used the saved file name if it exists
            if self.file_name_saved is None:
                raise ValueError("No file name specified to dump the json into!")
            else:
                # Use the saved file
                return self.file_name_saved
