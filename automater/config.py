import os
class Subconfig(object):
	pass


# System settings will be related to the computer or server that the quantum espresso code in running on.
system = Subconfig()
system.os = 'Linux' #For determining terminal commands
system.hasDisplay = 'DISPLAY' in os.environ # Setting to false will prevent plotting routines from trying to display figures
system.qeDir = r'/opt/qe-6.1' # The directory of the quantum espresso installation. 
system.pgmDir = r'/opt/qe-6.1/bin' # The bin directory of the quantum espresso installation
system.pseudoDir = r'/home/buck/qe/pseudo' # The directory where pseudopotential files are stored.
system.structureDir = r'/home/buck/qe/structures' # The directory where all structure files will be saved and loaded (by default) from 
system.projectDir = r'/home/buck/qe/projects' #The directory where all qe input and output files will be stored

# qeBot can email you when certain calculations have completed. These email settings are the settings of your chosen SMTP server. 
email = Subconfig()
email.able = True
email.host = 'dallas170.arvixeshared.com'
email.username = 'qebot@terratribe.org'
email.isSsl = True
email.port = 465
email.password = 'qebot238twO'
email.fromAddress = 'qebot@terratribe.org'
email.to = ['buck54321@gmail.com'] # List of email addresses to send reports

twilio = Subconfig()
twilio.sid = "AC98916e97429d270f3f08787f69ce9be7" # Your Account SID from www.twilio.com/console
twilio.token = "665c8bf60fb12d7ab49047e72d4bd1bd" # Your Auth Token from www.twilio.com/console
twilio.number = "+13124873268" # Your number from twilio. For U.S. numbers be sure to include the "+1" at the beginning of the number
twilio.to = "+18155662246" # The number to receive text reports. For U.S. numbers be sure to include the "+1" at the beginning of the number