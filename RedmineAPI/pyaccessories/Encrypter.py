import base64


class Encrypter(object):
    def __init__(self, path, name, key='Sixteen byte key'):
        """
        :param path: Path where you want to store/load the credential.
        :param name: Name of the credential for asking the user.
        :param key: 16 characters to use as a key for the encryption (can be left as default).
        """
        self.key = key
        self.path = path
        self.name = name

    def load(self):
        """
        :return: loaded unencrypted credential as string.
        """
        try:
            f = open(self.path, "r")
            return self.__decode(self.key, f.read())
        except FileNotFoundError:
            f = open(self.path, "w")
            cred = input("Please enter your %s: " % self.name)
            f.write(self.__encode(self.key, cred))
            return cred

    @staticmethod
    def __encode(key, string):
        encoded_chars = []
        for i in range(len(string)):
            key_c = key[i % len(key)]
            encoded_c = chr(ord(string[i]) + ord(key_c) % 256)
            encoded_chars.append(encoded_c)
        encoded_string = "".join(encoded_chars)
        encoded_string = bytes(encoded_string, "utf-8")

        return base64.urlsafe_b64encode(encoded_string).decode('utf-8')

    @staticmethod
    def __decode(key, string):
        decoded_chars = []
        string = base64.urlsafe_b64decode(string).decode('utf-8')
        for i in range(len(string)):
            key_c = key[i % len(key)]
            encoded_c = chr(abs(ord(str(string[i]))
                                - ord(key_c) % 256))
            decoded_chars.append(encoded_c)
        decoded_string = "".join(decoded_chars)

        return decoded_string
