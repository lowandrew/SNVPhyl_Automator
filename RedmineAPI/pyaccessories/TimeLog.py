import datetime as dt


class Timer(object):
    def __init__(self, log_file=None):
        self.start_time = dt.datetime.now()
        self.colour = None
        self.logger = None
        # If there is a logfile specified
        self.logging = type(log_file) == str

        if self.logging:
            import logging
            self.logger = logging.getLogger('scope.name')

            file_log_handler = logging.FileHandler(log_file)
            self.logger.addHandler(file_log_handler)

            # nice output format
            formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
            file_log_handler.setFormatter(formatter)

            self.logger.setLevel('DEBUG')

    def time_str(self, do_colour=True):
        if type(self.colour) is not str or not do_colour:
            return "[Elapsed time: %.2f] " % (dt.datetime.now() - self.start_time).total_seconds()
        else:
            return self.colour + "[Elapsed time: %.2f] " % (dt.datetime.now() - self.start_time).total_seconds() + \
                   "\033[0m"

    def reset(self):
        self.start_time = dt.datetime.now()

    def time_print(self, to_print):
        if self.logging:
            import logging
            self.logger.info(self.time_str(do_colour=False) + str(to_print))
        print(self.time_str() + str(to_print))

    def set_colour(self, colour=None):
        if 30 <= colour <= 37:
            self.colour = "\033[1;" + str(colour) + ";0m"
        elif colour is None:
            self.colour = None
        else:
            print("[Warning] Invalid colour!")
