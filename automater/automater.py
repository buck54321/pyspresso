from compatibility import *
import qeio
import pointgroup as pg
import sys
import pickle
import config
import os
import re

''' This section is for taking command line calls to kill or restart a calculation with the following syntax
	Kill:		automater.py -x PROJECT 
	Restart:	automater.py -r PROJECT 
	where the value of PROJECT corresponds to the calculation you want to kill.
'''
if len(sys.argv) > 1:
	if sys.argv[1] == '-x':
		if len(sys.argv) > 2:
			project = sys.argv[2]
			testArray = project.split('/')
			io = qeio.io(project)
			io.prefix = '-'.join(list(reversed(testArray)))+'-res'
			io.pause()
			exit()
		else:
			print('In order to end a calculation, you must provide a second command line argument corresponding to the project.')
			exit()
	elif sys.argv[1] == '-r':
		if len(sys.argv) > 2:
			filepath = os.path.join(config.system.qeDir, 'projects', sys.argv[2], 'io.pkl')
			if os.path.isfile(filepath):
				pickleFile = open(filepath, 'r')
				io = pickle.load(pickleFile)
				io.restartCalc()
				exit()
		else:
			print('In order to restart a calculation, you must provide a second command line argument corresponding to the project.')
			exit()
	else:
		print('Unkown command line argument detected. Automater will run as normal.')
		
# END COMMAND LINE START/STOP

#Set mpi parameters here [np,ni,nk,nb,nt,nd]
mpi=[32,1,4,1,2,1]
#mpi=False


'''
# CdTe primitive cell - PBE
io = qeio.io('CdTe-primitive', ecutwfc = 100, ecutrho = 800, kxyz = [10,10,10], mpi=mpi) # a = 4.701, ecutwfc = 100, ecutrho = 800, kxyz = [10,10,10]
io.importStructure('CdTe_mp-406_primitive.cif', format='cif')
io.setPseudoType('pz-mt_fhi')
io.subTest('test')
'''

'''
io.importStructure('CdTe_mp-406_primitive.cif', format='cif')
io.setA(4.701)
io.addElement('Cd',112.4,'Cd.pbe-mt_fhi.UPF')
io.addElement('Te', 127.6, 'Te.pbe-mt_fhi.UPF')
io.hubbardU = {'Cd':15.0}
'''
'''
# Cd primitive cell
#cellParams = [[1.,0.,0.],[-0.5,0.866,0.],[0.,0.,1.89]]
io = qeio.io('Cd-primitive', ecutwfc = 100, ecutrho = 800, kxyz = [24,24,24], mpi=mpi)
io.importStructure('Cd-primitive.cif', format='cif')
io.addElement('Cd',mass=112.4, pseudoFile='Cd.pw-mt_fhi.UPF')
#io.addAtom('Cd',[0.,0.,0.],[0,0,0])
#io.addAtom('Cd',[0.,0.5744,0.945],[1,1,1])
io.occupations = 'smearing'
io.degauss = 0.02
'''
'''
# Te primitive cell
io = qeio.io('Te-primitive', ecutwfc = 58, ecutrho = 464, kxyz = [6,6,6], mpi=mpi)
io.importStructure('Te_mp-19_primitive.cif', format='cif')
io.addElement('Te', 127.6, 'Te.pbe-mt_fhi.UPF')
io.occupations = 'smearing'
io.degauss = 0.02
'''
'''
#As primitive cell
io = qeio.io('As-primitive', ecutwfc = 44, ecutrho = 272, kxyz = [14,14,14], mpi=mpi)
io.importStructure('As_mp-11_primitive.cif', format='cif')
io.addElement('As',pseudoFile="As.pbesol-n-kjpaw_psl.0.2.UPF")
io.occupations = 'smearing'
io.degauss = 0.02
'''
'''
#Cd3As2 primitive cell
io = qeio.io('Cd3As2-primitive', ecutwfc = 100, ecutrho = 800, kxyz = [6,6,6], mpi=mpi)
io.importStructure('Cd3As2_mp-1372_primitive.cif', format='cif')
io.addElement('As',pseudoFile="As.pbe-mt_fhi.UPF")
io.addElement('Cd',pseudoFile='Cd.pbe-mt_fhi.UPF')
'''
'''
# CdTe:As conventional cell
io = qeio.io('CdTe-As-conventional', a = 6.5547, ecutwfc = 100, ecutrho = 800, kxyz = [4,6,6], mpi=mpi)
io.importStructure('pbe-mt_fhi-x211-CdTe-As-conventional-res-relax.cif', format='cif')
io.addElement('Cd',112.4,'Cd.pbe-mt_fhi.UPF')
io.addElement('Te', 127.6, 'Te.pbe-mt_fhi.UPF')
io.addElement('As', 74.9, 'As.pbe-mt_fhi.UPF')
io.charge = -1.0
#io.setA(6.5547)
#io.repeatStructure(0,1)
#io.replaceAtom('Te','As')
io.occupations = 'smearing'
io.degauss = 0.02
'''

io = qeio.io('CdTe-primitive', ecutwfc = 86, ecutrho = 688, kxyz = [4,4,4], mpi=mpi)
io.optimize.plot=True
io.subTest('PBE0')
io.importStructure('CdTe_mp-406_primitive.cif', format='cif')
io.setPseudoType('pbe-(d)*n-kjpaw_psl.0.2(.2)*')
io.input_dft = 'PBE0'
ecutwfc = io.optimize.ecutwfc(30, 200, 4, rhoRatio = 4)
if not ecutwfc:
	exit('ecutwfc not converged. Exiting convergence routine')
io.ecutwfc = ecutwfc
kGrid = io.optimize.kGrid([4,4,4], 20, stepsize = [2,2,2])
if not kGrid:
	exit('k-grid not converged. Exiting convergence routine.')
io.kxyz = kGrid
print('k-grid converged at [%s]' % ','.join([str(x) for x in kGrid]))
io.optimize.relax(20, vc=False)
io.loadRelaxedStructure()
io.setA(io.optimize.latticeConstant(io.crystal.a*0.9, io.crystal.a*1.1, 10))
io.optimize.bandPlot([pg.symPts.fcc.L, pg.symPts.fcc.gamma, pg.symPts.fcc.X], 50, bandLabels = ['L', '$\Gamma$', 'X'])
io.optimize.dosPlot(Ecenter=False, Espread=False, Emin=-10, Emax=9, deltaE = 0.02)
io.rmDir('results')
exit()
'''
pseudoTypes = ['pbe-mt_fhi', 'pbe-(d)*n-kjpaw_psl.0.2(.2)*', 'pbe-(d)*n-rrkjus_psl.0.2(.2)*', 'pbesol-(d)*n-kjpaw_psl.0.2(.2)*', 'pbesol-(d)*n-rrkjus_psl.0.2(.2)*', 'pz-mt_fhi','pz-(d)*n-kjpaw_psl.0.2(.2)*', 'pz-(d)*n-rrkjus_psl.0.2(.2)*'] #, 'rel-pbe-(d)*n-kjpaw_psl.0.2(.2)*', 'rel-pbe-(d)*n-rrkjus_psl.0.2(.2)*', 'rel-pbesol-(d)*n-kjpaw_psl.0.2(.2)*', 'rel-pbesol-(d)*n-rrkjus_psl.0.2(.2)*', 'rel-pz-(d)*n-kjpaw_psl.0.2(.2)*', 'rel-pz-(d)*n-rrkjus_psl.0.2(.2)*']
#pseudoTypes = ['rel-pbe-(d)*n-kjpaw_psl.0.2(.2)*', 'rel-pbe-(d)*n-rrkjus_psl.0.2(.2)*', 'rel-pbesol-(d)*n-kjpaw_psl.0.2(.2)*', 'rel-pbesol-(d)*n-rrkjus_psl.0.2(.2)*', 'rel-pz-(d)*n-kjpaw_psl.0.2(.2)*', 'rel-pz-(d)*n-rrkjus_psl.0.2(.2)*']
hubbardUs = [False, 5., 10., 15.]
for ptype in pseudoTypes:
	for U in hubbardUs:
		io = qeio.io('CdTe-primitive', ecutwfc = 100, ecutrho = 800, kxyz = [6,6,6], mpi=mpi) # a = 4.701, ecutwfc = 100, ecutrho = 800, kxyz = [10,10,10]
		io.waitAtImages = False
		io.importStructure('CdTe_mp-406_primitive.cif', format='cif')
		if re.search('^rel', ptype):
			io.makeRelativistic()
		typeString = re.sub('[^a-zA-Z0-9_-]', '', ptype)
		io.subTest(typeString)	
		io.setPseudoType(ptype)
		if U:
			io.hubbardU = {'Cd':U}
			io.subTest('hubbardU%i' % int(U), True)
		else:
			io.hubbardU = False
			io.subTest('noHubbard', True)
		print('------------------------------------------------')
		print('Starting calculations for pseudopotential %s' % typeString)
		if U:
			print('Using hubbard U correction on Cd of %.1f eV' % U)
		io.importStructure(io.postfix+'-relax.cif')
		try: 
			dictionary = io.readDictionary()
			io.kxyz = [int(x) for x in dictionary['k-grid'].split('-')]
			io.ecutwfc = int(dictionary['ecutwfc'])
			io.ecutrho = int(dictionary['ecutrho'])
			io.setA(float(dictionary['lattice-constant']))
		except:
			io.log(io.highlight('Error parsing dictionary for %s' % io.postfix))
			continue
		io.optimize.bandPlot([pg.symPts.fcc.L, pg.symPts.fcc.gamma, pg.symPts.fcc.X], 50, bandLabels = ['L', '$\Gamma$', 'X'])


		ecutwfc = io.optimize.ecutwfc(30, 200, 4, rhoRatio = 8)
		if not ecutwfc:
			exit('ecutwfc not converged. Exiting convergence routine')
		io.ecutwfc = ecutwfc
		kGrid = io.optimize.kGrid([4,4,4], 20, stepsize = [2,2,2])
		if not kGrid:
			exit('k-grid not converged. Exiting convergence routine.')
		io.kxyz = kGrid
		print('k-grid converged at [%s]' % ','.join([str(x) for x in kGrid]))
		io.optimize.relax(20, vc=False)
		io.loadRelaxedStructure()
		io.setA(io.optimize.latticeConstant(io.crystal.a*0.9, io.crystal.a*1.1, 10))
		io.optimize.relax(20, vc=True)
		io.loadRelaxedStructure()
		io.optimize.bandPlot([pg.symPts.fcc.L, pg.symPts.fcc.gamma, pg.symPts.fcc.X], 50, bandLabels = ['L', '$\Gamma$', 'X'])
		io.optimize.dosPlot(Ecenter=False, Espread=False, Emin=-10, Emax=9, deltaE = 0.02)
		io.rmDir('results')
'''


#io.optimize.ecutwfc(80, 200, 4, rhoRatio = 8)
#io.optimize.ecutrho(600, 1000, 16)
#io.optimize.kGrid([4,4,4], 10, stepsize = [2,2,2])

#io.optimize.latticeConstant(io.crystal.a*0.9, io.crystal.a*1.1, 10)

#io.optimize.relax(20, vc=False, constraints = False, cellDoFree=False)

#io.optimize.quickRun()

#io.optimize.bandPlot([pg.symPts.fcc.L, pg.symPts.fcc.gamma, pg.symPts.fcc.X], 50, bandLabels = ['L', '$\Gamma$', 'X'], nbnd=18)
#io.plotBands()

#io.optimize.dosPlot(Ecenter=False, Espread=False, Emin=-10, Emax=9, deltaE = 0.01, nbnd=18)

'''
io.optimize.relax(20, vc=False)
io.loadRelaxedStructure()
io.optimize.relax(20, vc=True)
'''

#io.exportData()





'''
PSEUDOPOTENTIALS
Si.blyp-hgh.UPF
Si.pz-vbc.UPF
Si.pw91-n-van.UPF

As.pbe-mt_fhi.UPF
As.pbe-n-kjpaw_psl.0.2.upf
As.pbe-n-kjpaw_psl.0.2.UPF
As.pbe-n-rrkjus_psl.0.2.UPF
As.pbe-n-van.UPF
As.pbesol-n-kjpaw_psl.0.2.UPF
As.pbesol-n-rrkjus_psl.0.2.UPF
As.pw91-n-van.UPF
As.pz-bhs.UPF
As.pz-mt_fhi.UPF
As.pz-n-kjpaw_psl.0.2.UPF
As.pz-n-rrkjus_psl.0.2.UPF
As.rel-pbe-n-kjpaw_psl.0.2.UPF
As.rel-pbe-n-rrkjus_psl.0.2.UPF
As.rel-pbesol-n-kjpaw_psl.0.2.UPF
As.rel-pbesol-n-rrkjus_psl.0.2.UPF
As.rel-pz-n-kjpaw_psl.0.2.UPF
As.rel-pz-n-rrkjus_psl.0.2.UPF

Cd.pbe-dn-kjpaw_psl.0.2.UPF
Cd.pbe-dn-rrkjus_psl.0.2.UPF
Cd.pbe-mt_EM.UPF
Cd.pbe-mt_fhi.UPF
Cd.pbe-n-van.UPF
Cd.pbesol-dn-kjpaw_psl.0.2.UPF
Cd.pbesol-dn-rrkjus_psl.0.2.UPF
Cd.pw91-n-van.UPF
Cd.pw-mt_fhi.UPF
Cd.pz-dn-kjpaw_psl.0.2.UPF
Cd.pz-dn-rrkjus_psl.0.2.UPF
Cd.pz-mt_fhi.UPF
Cd.rel-pbe-dn-kjpaw_psl.0.2.UPF
Cd.rel-pbe-dn-rrkjus_psl.0.2.UPF
Cd.rel-pbesol-dn-kjpaw_psl.0.2.UPF
Cd.rel-pbesol-dn-rrkjus_psl.0.2.UPF
Cd.rel-pz-dn-kjpaw_psl.0.2.UPF
Cd.rel-pz-dn-rrkjus_psl.0.2.UPFTe.pw-mt_fhi.UPF

Te.pbe-dn-kjpaw_psl.0.2.2.UPF
Te.pbe-dn-rrkjus_psl.0.2.2.UPF
Te.pbe-mt_fhi.UPF
Te.pbe-rrkj.UPF
Te.pbesol-dn-kjpaw_psl.0.2.2.UPF
Te.pbesol-dn-rrkjus_psl.0.2.2.UPF
Te.pw-mt_fhi.UPF
Te.pz-bhs.UPF
Te.pz-dn-kjpaw_psl.0.2.2.UPF
Te.pz-dn-rrkjus_psl.0.2.2.UPF
Te.pz-mt_fhi.UPF
Te.rel-pbe-dn-kjpaw_psl.0.2.2.UPF
Te.rel-pbe-dn-rrkjus_psl.0.2.2.UPF
Te.rel-pbesol-dn-kjpaw_psl.0.2.2.UPF
Te.rel-pbesol-dn-rrkjus_psl.0.2.2.UPF
Te.rel-pz-dn-kjpaw_psl.0.2.2.UPF
Te.rel-pz-dn-rrkjus_psl.0.2.2.UPF
'''
'''
# Si primitive cell
io = qeio.io('si-primitive', a = 5.468, ecutwfc = 28, ecutrho = 112, kxyz = [6,6,6], ibrav=2, mpi=mpi)
io.addElement('Si',28.085,'Si.pw91-n-van.UPF')
io.addAtom('Si',[0.0,0.0,0.0])
io.addAtom('Si',[0.25,0.25,0.25])
'''
'''
# Si conventional cell
cellParams = [[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]]
io = qeio.io('si-111', a = 5.466653, ecutwfc = 26, ecutrho = 112, kxyz = [14,14,14], ibrav=0, cellParams=cellParams, mpi=mpi)
io.addElement('Si',28.085,'Si.pw91-n-van.UPF')
vectors =  [[0.0,0.0,0.0],[0.5,0.5,0.0],[0.5,0.0,0.5],[0.0,0.5,0.5]]
for atom in vectors:
	io.addAtom('Si',atom)
#for atom in vectors[1:]:
for atom in vectors:
	io.addAtom('Si',[x+0.25 for x in atom])
'''
'''
# 2x1x1 Si supercell
cellParams = False
cellParams = [[2.0,0.0,0.0], [0.0,1.0,0.0], [0.0,0.0,1.0]]
io = qeio.io('si-211-vacancy', a = 5.466653, ecutwfc = 32, ecutrho = 140, kxyz = [8,10,10], ibrav=0, mpi=mpi, cellParams = cellParams)
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
io = qeio.io('si-221', a = 5.466653, ecutwfc = 34, ecutrho = 158, kxyz = [7,7,7], ibrav=0, mpi=mpi, cellParams = cellParams)
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
io = qeio.io('si-222', a = 5.466653, ecutwfc = 34, ecutrho = 158, kxyz = [7,7,7], ibrav=0, mpi=mpi, cellParams = cellParams)
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
io = qeio.io('inas-primitive', a = 6.1914, ecutwfc = 25, ecutrho = 100, kxyz = [4,4,4], ibrav=2)
io.addElement('In', 114.818, 'In.pbe-dn-kjpaw_psl.0.2.2.upf')
io.addElement('As', 74.922, 'As.pbe-n-kjpaw_psl.0.2.upf')
io.addAtom('In',[0.0,0.0,0.0])
io.addAtom('As',[0.25,0.25,0.25])
'''
'''
# CdTe conventional cell
io = qeio.io('CdTe-conventional', a = 6.48, ecutwfc = 48, ecutrho = 196, kxyz = [4,4,4], mpi=mpi)
io.importStructure('CdTe_mp-406_conventional_standard.cif', format='cif')
io.addElement('Cd',112.4,'Cd.pbesol-dn-kjpaw_psl.0.2.UPF')
io.addElement('Te', 127.6, 'Te.pbesol-dn-kjpaw_psl.0.2.2.UPF')
io.printStructure()
'''
