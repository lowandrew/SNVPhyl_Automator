# pyaccessories
All my python accessory classes

## Installation

```console
git clone https://github.com/devonpmack/pyaccessories.git
cd dist
pip install py-accessories-*.tar
```

### Submodule
```console
git submodule add https://github.com/devonpmack/pyaccessories.git
```

#### Encrypter.py
Instantiate the encrypter with `var = Encrypter(path, name)`
- path: Path where you want to store/load the credential.
- Name of the credential for asking the user. eg. Api Key
- key: (optional) 16 characters to use as a key for the encryption (by default will use 'Sixteen byte key').

Then use the `Encrypter.load()` method to attempt to load your encrypted credential from a file. If the file doesn't exist the program will ask for the credential and save it to the file encrypted.
##### Usage
```python
from pyaccessories.Encrypter import Encrypter
e = Encrypter('pass.txt', 'Api Key')
print(e.load())
```

#### SaveLoad.py
```python
from pyaccessories.SaveLoad import SaveLoad

# Create a new SaveLoad instance which loads the file config.json, creating it
# if it doesn't already exist
loader = SaveLoad('config.json', create=True)

# Loads the value 'output' from the config file and stores it in the variable output_file
# if output isn't in the config file then ask the user to input it. To turn this ask off use
# ask=False and then it will throw an error if the config option isn't there
output_file = loader.get('output', default='/home/out.txt')

# Use the loaded config option and write to the file
open(output_file, 'w').write('Hello World')
```

#### TimeLog.py
```python
from pyaccessories.TimeLog import Timer

# Create a new instance of Timer which will log all time_prints to log.txt
# If no log file is specified it will not log
t = Timer(log_file='log.txt')

# Set the colour to lime green
t.set_colour(30)

# Print the message hi
t.time_print('Hi')
```
