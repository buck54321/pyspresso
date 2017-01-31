import numpy as np
import matplotlib.pyplot as plt

dt = 0.01 # time step in seconds
length = 10.16 # length of plate in centimeters
width = 10.16 # width of plate in centimeters
height = 0.9525 # height of plate in centimeters
A = length*width #area of plate in square centimeters
Am = A*0.0001 # area of plate in square meteres
V = A*height # volume of plate in cubic centimeters
c = 0.900 # heat capacity of aluminum in J/gC
p = 2.7 # density of aluminum in g/cm^3
m = V*p # mass of the plate in grams
C = c*m # total heat capacity of the plate
W = 250.0 # power output of peltier coolers in watts
sigma = 5.67e-8 # steffan boltmann blackbody factor
maxDT = 50.0 # maximum temperature difference of peltier coolers
eAir = 10 # Air convection power loss coefficient
Tenv = 300 # environmental air temperature
Tplate = 300 # initial temperature of the cold plate
TplateNoLoss = 300
TplateAirOnly = 300
TplateBBonly = 300

def e(dT):
	"""For a given temperature difference, caluclate the efficiency of the peltier coolers"""
	if dT > 50:
		dt = 50.0
	return ((50.0 - dT)/50.0)*0.02

def Tin(Th, Tc, W=W):
	"""The total temperature per time from input power"""
	return e(Th-Tc)*W/C

def Tbb(Th, Tc):
	"""Total temperature change per time from blackbody radiation"""
	return sigma*Am*(Th**4 - Tc**4)/C

def Tair(Th, Tc):
	return eAir*Am*(Th - Tc)/C

t = 0
tlog = []
Tlog = []
TnoLossLog = []
TairOnlyLog = []
TbbOnlyLog = []
i = 0
while True:
	TplateNoLoss += -1*Tin(Tenv,TplateNoLoss)*dt
	TplateAirOnly += -1*Tin(Tenv,TplateAirOnly)*dt + Tair(Tenv,TplateAirOnly)*dt
	TplateBBonly += -1*Tin(Tenv,TplateBBonly)*dt + Tbb(Tenv,TplateAirOnly)*dt
	Tplate += -1*Tin(Tenv,Tplate)*dt + Tbb(Tenv,Tplate)*dt +  Tair(Tenv,Tplate)*dt
	#print 'Tplate: %f' % Tplate
	t += dt
	if i%10000 is 0:
		tlog.append(t/60.0)
		Tlog.append(Tplate - 273)
		TnoLossLog.append(TplateNoLoss - 273)
		TairOnlyLog.append(TplateAirOnly - 273)
		TbbOnlyLog.append(TplateBBonly - 273)
	i += 1
	if i > 500000:
		break

plt.ion()
fig = plt.figure()
fig.patch.set_facecolor('white')
#fig.xlabel = 'time (min)'
#fig.ylabel = 'tempearature (C)'
ax = plt.gca()
ax.ticklabel_format(useOffset=False)
ax.set_ylabel('Temperature (C)', fontsize=18)
ax.set_xlabel('Time (min)', fontsize=18)
Tmarkers = plt.plot(tlog,Tlog)
plt.setp(Tmarkers, color='k', marker='s', linestyle = 'none', markersize = 6.0, label = 'Air and BB')
TnoLossMarkers = plt.plot(tlog,TnoLossLog)
plt.setp(TnoLossMarkers , color='r', marker='o', linestyle = 'none', markersize = 6.0, label = 'No losses')
TairOnlyMarkers = plt.plot(tlog,TairOnlyLog)
plt.setp(TairOnlyMarkers , color='b', marker='v', linestyle = 'none', markersize = 6.0, label = 'Air Only')
TbbOnlyMarkers = plt.plot(tlog,TbbOnlyLog)
plt.setp(TbbOnlyMarkers , color='g', marker='+', linestyle = 'none', markersize = 6.0, label = 'BB Only')
plt.suptitle('Peltier Cold Plate Temperature Simulation', fontsize=20)
plt.legend(loc="upper right", fancybox=True, ncol=1, fontsize=13, numpoints = 1)
plt.show()
raw_input('Press enter to continue')




