import sys
v = sys.version_info[0]

if v >= 3:
	def raw_input(str):
		return input(str)

def iteritems(d):
	if v >=3:
		return d.items()
	return d.iteritems()

