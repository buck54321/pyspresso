import config
import smtplib
from os.path import basename
from email.mime.application import MIMEApplication
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.utils import COMMASPACE, formatdate

from twilio import TwilioRestException
from twilio.rest import TwilioRestClient

def email(subject, message, to=config.email.to, files=False):
	if config.email.isSsl:
		emailer = smtplib.SMTP_SSL(config.email.host, config.email.port)
	else:
		emailer = smtplib.SMTP(config.email.host, config.email.port)
		emailer.ehlo()
		emailer.starttls()
	emailer.login(config.email.username, config.email.password)
	msg = MIMEMultipart()
	msg['From'] = config.email.fromAddress
	msg['To'] = COMMASPACE.join(to)
	msg['Date'] = formatdate(localtime=True)
	msg['Subject'] = subject
	msg.attach(MIMEText(message))
	for f in files or []:
		with open(f, "rb") as fil:
			fileName = basename(f)
			part = MIMEApplication(fil.read(), Name=fileName)
			part['Content-Disposition'] = 'attachment; filename="%s"' % fileName
			msg.attach(part)
	emailer.sendmail(config.email.fromAddress, config.email.to, msg.as_string())
	emailer.quit()
	
	

def text(msg, to=config.twilio.to):
	client = TwilioRestClient(config.twilio.sid, config.twilio.token)
	try:
		message = client.messages.create(body = msg,
			to=to,    # Replace with your phone number
			from_= config.twilio.number) # Replace with your Twilio number
		return 'Great'
	except TwilioRestException as e:
		return e