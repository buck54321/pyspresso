import qeio
import optimize
import sys
import pickle
import config
import os

''' This section is for taking command line calls to kill or restart a calculation with the following syntax
	Kill:		automater.py -x PROJECT 
	Restart:	automater.py -r PROJECT 
	where the value of PROJECT corresponds to the calculation you want to kill.
'''
if len(sys.argv) > 1:
	if sys.argv[1] == '-x':
		if len(sys.argv) > 2:
			project = sys.argv[2]
			io = qeio.io(project)
			io.pause()
			exit()
		else:
			print 'In order to end a calculation, you must provide a second command line argument corresponding to the project.'
			exit()
	elif sys.argv[1] == '-r':
		if len(sys.argv) > 2:
			filepath = os.path.join(config.system['qeDir'], 'projects', sys.argv[2], 'io.pkl')
			if os.path.isfile(filepath):
				pickleFile = open(filepath, 'r')
				io = pickle.load(pickleFile)
				io.restartCalc()
				exit()
		else:
			print 'In order to restart a calculation, you must provide a second command line argument corresponding to the project.'
			exit()
	else:
		print 'Unkown command line argument detected. Automater will run as normal.'
		
# END COMMAND LINE START/STOP

#Set mpi parameters here [np,ni,nk,nb,nt,nd]
#mpi=[16,1,4,1,2,1]
mpi=False

'''
# Si primitive cell
io = qeio.io('si-primitive', a = 5.466653, eCutWfc = 26, eCutRho = 112, kxyz = [14,14,14], ibrav=2)
io.addElement('Si',28.085,'Si.pw91-n-van.UPF')
io.addAtom('Si',[0.0,0.0,0.0])
io.addAtom('Si',[0.25,0.25,0.25])
'''

'''
# Si conventional cell
io = qeio.io('si-111', a = 5.466653, eCutWfc = 26, eCutRho = 112, kxyz = [14,14,14], ibrav=1)
io.addElement('Si',28.085,'Si.pw91-n-van.UPF')
vectors = [[0.0,0.0,0.0],[0.5,0.5,0.0],[0.5,0.0,0.5],[0.0,0.5,0.5]]
for atom in vectors:
	io.addAtom('Si',atom)
for atom in vectors[1:]:
#for atom in vectors:
	io.addAtom('Si',[x+0.25 for x in atom])
'''

'''
# 2x1x1 Si supercell
cellParams = False
cellParams = [[2.0,0.0,0.0], [0.0,1.0,0.0], [0.0,0.0,1.0]]
io = qeio.io('si-211-vacancy', a = 5.466653, eCutWfc = 32, eCutRho = 140, kxyz = [8,10,10], ibrav=0, mpi=mpi, cellParams = cellParams)
io.addElement('Si',28.085,'Si.pw91-n-van.UPF')
vectors = [[0.0,0.0,0.0],[0.5,0.5,0.0],[0.5,0.0,0.5],[0.0,0.5,0.5]]
scTranslations = [[1.0,0.0,0.0]]
for atom in vectors:
	io.addAtom('Si',atom)
for atom in vectors[1:]:
#for atom in vectors:
	io.addAtom('Si',[x+0.25 for x in atom])
for translation in scTranslations:
	for atom in vectors:
		io.addAtom('Si',[x+y for x,y in zip(atom,translation)])
		io.addAtom('Si',[x+y+0.25 for x,y in zip(atom,translation)])
'''

'''
# 2x2x1 Si supercell
cellParams = False
cellParams = [[2.0,0.0,0.0], [0.0,2.0,0.0], [0.0,0.0,1.0]]
io = qeio.io('si-221', a = 5.466653, eCutWfc = 34, eCutRho = 158, kxyz = [7,7,7], ibrav=0, mpi=mpi, cellParams = cellParams)
io.addElement('Si',28.085,'Si.pw91-n-van.UPF')
vectors = [[0.0,0.0,0.0],[0.5,0.5,0.0],[0.5,0.0,0.5],[0.0,0.5,0.5]]
scTranslations = [[1.0,0.0,0.0], [0.0,1.0,0.0], [1.0,1.0,0.0]]
for atom in vectors:
	io.addAtom('Si',atom)
#for atom in vectors[1:]:
for atom in vectors:
	io.addAtom('Si',[x+0.25 for x in atom])
for translation in scTranslations:
	for atom in vectors:
		io.addAtom('Si',[x+y for x,y in zip(atom,translation)])
		io.addAtom('Si',[x+y+0.25 for x,y in zip(atom,translation)])
'''

'''
# 2x2x2 Si supercell
cellParams = False
cellParams = [[2.0,0.0,0.0], [0.0,2.0,0.0], [0.0,0.0,2.0]]
io = qeio.io('si-222', a = 5.466653, eCutWfc = 34, eCutRho = 158, kxyz = [7,7,7], ibrav=0, mpi=mpi, cellParams = cellParams)
io.addElement('Si',28.085,'Si.pw91-n-van.UPF')
vectors = [[0.0,0.0,0.0],[0.5,0.5,0.0],[0.5,0.0,0.5],[0.0,0.5,0.5]]
scTranslations = [[1.0,0.0,0.0], [0.0,1.0,0.0], [1.0,1.0,0.0], [0.0,0.0,1.0], [1.0,0.0,1.0], [0.0,1.0,1.0], [1.0,1.0,1.0]]
for atom in vectors:
	io.addAtom('Si',atom)
#for atom in vectors[1:]:
for atom in vectors:
	io.addAtom('Si',[x+0.25 for x in atom])
for translation in scTranslations:
	for atom in vectors:
		io.addAtom('Si',[x+y for x,y in zip(atom,translation)])
		io.addAtom('Si',[x+y+0.25 for x,y in zip(atom,translation)])
'''

'''
# Indium Arsenide primitive cell
io = qeio.io('inas-primitive', a = 6.1914, eCutWfc = 25, eCutRho = 100, kxyz = [4,4,4], ibrav=2)
io.addElement('In', 114.818, 'In.pbe-dn-kjpaw_psl.0.2.2.upf')
io.addElement('As', 74.922, 'As.pbe-n-kjpaw_psl.0.2.upf')
io.addAtom('In',[0.0,0.0,0.0])
io.addAtom('As',[0.25,0.25,0.25])
'''


# CdTe primitive cell
io = qeio.io('cdte-primitive', a = 6.48, eCutWfc = 49, eCutRho = 186, kxyz = [12,12,12], ibrav=2)
io.addElement('Cd',112.4,'Cd.pbe-dn-rrkjus_psl.0.2.UPF')
io.addElement('Te', 127.6, 'Te.pbe-dn-rrkjus_psl.0.2.2.UPF')
io.addAtom('Cd',[0.0,0.0,0.0])
io.addAtom('Te',[0.25,0.25,0.25])


class Params(object):
	pass
optParams = Params()
optParams.startingwfc = False # Setting to 'file' will use wavefunctions from previous calculation as initial wavefunctions in new calculation. 
optParams.startingpot = False # Setting to 'file' will use charge density from previous calculation as initial charge density for new calculation.
optParams.plot = True # Setting to False will prevent creation of plots and plot files. 
optParams.replot = False # Setting to True will prevent recalculation, and simply parse data from previous 
optParams.restart = False # Setting to True will put the calculation in restart mode, which picks up from a previously paused calculation. See QE documentation on how to pause a calculation. 
optParams.email = False # Setting to True will cause an email to be sent at the end of the optimization. Email (SMTP settings) must be configured in configure.py file
optParams.text = False # Setting to True will cause an text message to be sent at the end of the optimization. Text message (Twilio SMS settings) must be configured in configure.py file
optParams.convEV = 0.0001 # Convergence threshold. If three consecutive iterations of an optimization differ by less convEV, the optimization will halt. 

#optimize.quickRun(io, optParams)
optimize.latticeConstant(io, 6.3, 6.7, 10, optParams)
#optimize.kGrid(io, [8,8,8], 10, optParams, stepsize = [2,2,2])
#optimize.cutWfc(io, 53, 61, 2, optParams, rhoRatio = 4)
#optimize.cutRho(io, 180, 232, 2, optParams)
#optimize.relax(io, 3, optParams, vc=False)
io.exportData()



'''
I can use this section to store some shit
Si.pz-vbc.UPF
Si.pw91-n-van.UPF
Cd.pbe-dn-rrkjus_psl.0.2.UPF
Cd.pbesol-dn-rrkjus_psl.0.2.UPF
Te.pbe-dn-rrkjus_psl.0.2.2.UPF
Te.pbesol-dn-rrkjus_psl.0.2.2.UPF
'''
