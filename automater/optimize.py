from compatibility import *
import os
import numpy as np
import time
import copy
import config
from collections import OrderedDict


ry2ev = 13.605693
bohr2ang = 0.529177249

class Optimizer:
	"""Designed to be a class of subfunctions of Qeio for running crystal optimization scripts"""
	def __init__(self, qeio):
		self.qeio = qeio
		self.salvageNumber = 0
		self.startingwfc = False # Setting to 'file' will use wavefunctions from previous calculation as initial wavefunctions in new calculation. 
		self.startingpot = False # Setting to 'file' will use charge density from previous calculation as initial charge density for new calculation.
		self.plot = False # Setting to False will prevent creation of plots and plot files. 
		self.replot = False # Setting to True will prevent recalculation, and simply parse data from previous 
		self.restart = False # Setting to True will put the calculation in restart mode, which picks up from a previously paused calculation. See QE documentation on how to pause a calculation. 
		self.email = False # Setting to True will cause an email to be sent at the end of the optimization. Email (SMTP settings) must be configured in configure.py file
		self.text = False # Setting to True will cause an text message to be sent at the end of the optimization. Text message (Twilio SMS settings) must be configured in configure.py file
		self.convEV = 'auto' 

	def latticeConstant(self, low, high, steps):
		oglow, oghigh, ogsteps = (low, high, steps)
		self.qeio.calculation = 'scf'
		calcType = 'lattice'
		if self.convEV == 'auto':
			self.convEV == self.qeio.nat*0.001
		"""optimize lattice constant by """
		subDir = self.qeio.latticeFileDir
		opDir = os.path.join(self.qeio.workingDir, subDir)
		if not self.replot:
			if self.restart: 
				if not os.path.isdir(opDir):
					self.qeio.log('You are attempting to restart a calculation, but the previous results were not found. Set restart=False to start a new calculation.')
					return
				self.qeio.restart = True
			else:
				self.qeio.clearDir(subDir)
			latticeConstants = np.linspace(low, high, steps)
			self.qeio.log('Lattice calculation commenced %s (server time)' % time.asctime())
			self.qeio.startingwfc = self.startingwfc
			self.qeio.startingpot = self.startingpot
			for a in latticeConstants:
				low = a #Necessary in case optimization must be restarted
				steps += -1 #Necessary in case optimization must be restarted
				self.qeio.optimization = 'lattice'
				self.qeio.optLocals = locals()
				self.qeio.setA(a)
				self.qeio.write(calcType)
				opFilename = os.path.join(opDir, '%s-lattice-%.2f.out' % (self.qeio.postfix, a))#opDir+'/lattice-%f' % a
				result = self.qeio.runPw(opFilename, calcType)
				if result.errorEncountered:
					if result.errorType == 'loneVector':
						if self.salvageNumber < 5:
							self.salvageNumber += 1
							self.qeio.ecutrho = int(np.floor(float(self.qeio.ecutrho)*1.1))
							self.qeio.log(self.qeio.highlight('Lone vector error encountered. Restarting calculation with increased ecutrho = %i' % self.qeio.ecutrho))
							return self.latticeConstant(oglow, oghigh, ogsteps)
						else:
							self.qeio.log(self.qeio.highlight('Could not resolve lone vector error. Exiting lattice routine.'))
				#self.qeio.startingpot = 'file'
				#self.qeio.startingwfc = 'file'
				data = self.qeio.parseResults(subDir, log='lattice')
				data = sorted(data, key=lambda k: float(k['latticeConstant']))
				pt = data[-1]
				self.qeio.log('E(%s) = %s' % (pt['latticeConstant'], pt['energy']))
				self.qeio.log('Iteration wall time: %is' % pt['wallTime'])
				self.qeio.log('Date/Time: %s' % data[-1]['dateTime'])
			self.qeio.parseResults(subDir, log='lattice')
		if os.path.isdir(opDir):
			returnDict = self.qeio.fitLattice(plot=self.plot)
			if returnDict:
				self.qeio.log("Quadratic data fit yields a lattice constant of %f" % returnDict['quadFit'])
				self.qeio.log("Birch-Murnaghan EOS fit yields a lattice constant of %f" % returnDict['BMfit'])
				self.qeio.log("BM equation bulk modulus fit yields a bulk modulus of %f GPa" % returnDict['bulkMod'])
				self.qeio.logValue({'bulk-modulus(GPa)': returnDict['bulkMod'], 'lattice-constant(Angstrom)': returnDict['BMfit']})
				if self.plot and config.system.hasDisplay and self.qeio.waitAtImages:
					raw_input('Press any key to continue')
			else:
				self.qeio.log('Lattice fitting failed with error: %s' % self.qeio.lastError)
			if self.email or self.text:
				message = '%s calculation complete on project: %s' % ('lattice', self.qeio.project)
				files = returnDict['files'] or False
				self.qeio.report(message, text=text, email=email, files=files)
		else:
			self.qeio.log("Lattice files directory not found. Is replot set to True without first running calculations with replot=False?")
		data = self.qeio.parseResults(subDir)
		return returnDict['BMfit']
			
	def ecutwfc(self, low, high, stepsize, rhoRatio = 4):
		args = locals()
		args.pop('self',None)
		args['qeParam'] = 'ecutwfc'
		return self.cutoffs(**args)

	def ecutrho(self, low, high, stepsize):
		args = locals()
		args.pop('self',None)
		args['qeParam'] = 'ecutrho'
		return self.cutoffs(**args)

	def cutoffs(self, low, high, stepsize, qeParam, rhoRatio = 4):
		self.qeio.calculation = 'scf'
		if self.convEV == 'auto':
			self.convEV == self.qeio.nat*0.001
		subDir = self.qeio.ecutwfcFileDir if qeParam == 'ecutwfc' else self.qeio.ecutrhoFileDir
		opDir =  os.path.join(self.qeio.workingDir, subDir)
		if self.replot:
			if os.path.isdir(opDir):
				returnDict = self.qeio.fitCutoff(qeParam, subDir)
				message = '%s calculation complete on project: %s' % (qeParam, self.qeio.project)
				self.qeio.report(message, text=self.text, email=self.email, files=returnDict['files'])
				return
		if self.restart: 
			if not os.path.isdir(opDir):
				self.qeio.log('You are attempting to restart a calculation, but the previous results were not found. Set restart=False to start a new calculation.')
				return
			#self.qeio.restart = True
		else:
			self.qeio.clearDir(subDir)
		self.qeio.startingwfc = self.startingwfc
		self.qeio.startingpot = self.startingpot
		testE = low
		trigger = False
		triggered = False
		converged = False
		if self.convEV == 'auto':
			self.convEV = self.qeio.nat * 0.001
		self.qeio.log('%s calculation commenced %s (server time)' % (qeParam, time.asctime()))
		while testE <= high:
			self.qeio.optimization = qeParam
			self.qeio.optLocals = locals()
			opFilename = os.path.join(opDir, '%s-%s-%i.out' % (self.qeio.postfix, qeParam, testE))#opDir+'/%s-%f' % (qeParam, testE)
			if qeParam == 'ecutwfc':
				self.qeio.ecutwfc = testE
				if rhoRatio:
					self.qeio.ecutrho = rhoRatio*testE
			if qeParam == 'ecutrho':
				self.qeio.ecutrho = testE
			self.qeio.write(qeParam)
			self.qeio.runPw(opFilename, qeParam)
			#self.qeio.startingwfc = 'file'
			#self.qeio.startingpot = 'file'
			data = self.qeio.parseResults(subDir)
			#print([str(i[qeParam]) for i in data])
			data = sorted(data, key=lambda k: float(k[qeParam]))
			self.qeio.log('Total energy for %s=%i: %f' % (qeParam, testE, data[-1]['energy']))
			self.qeio.log('Total time: %is' % data[-1]['wallTime'])
			self.qeio.log('Date/time: %s' % data[-1]['dateTime'])
			if testE > low:
				deltaEV = data[-1]['energy'] - data[-2]['energy']
				self.qeio.log(self.qeio.highlight('%f eV change' % deltaEV))
				if abs(deltaEV) < self.convEV:
					if triggered:
						if qeParam == 'ecutwfc':
							self.qeio.ecutwfc = trigger
							self.qeio.ecutrho = trigger*rhoRatio
						if qeParam == 'ecutrho':
							self.qeio.ecutrho = trigger
						self.qeio.logValue({'ecutwfc(Ry)':self.qeio.ecutwfc, 'ecutrho(Ry)':self.qeio.ecutrho})
						self.qeio.log('%s converged at %i' % (qeParam, testE))
						converged = True
						break
					else:
						trigger = copy.copy(testE)
						triggered = True
			testE += stepsize
		self.qeio.parseResults(subDir, log=qeParam)
		if not converged:
			self.qeio.log('%s did not converge' % qeParam)
		returnDict = self.qeio.fitCutoff(qeParam, subDir, plot=self.plot)
		if self.email or self.text:
			if not returnDict:
				message = '%s calculation on %s project completed after fewer than 5 iterations' % (qeParam, self.qeio.project)
				self.qeio.report(message, text=self.text, email=self.email)
			else:
				message = '%s calculation complete on project: %s' % (qeParam, self.qeio.project)
				self.qeio.report(message, text=self.text, email=self.email, files=returnDict['files'])
		if converged:
			return trigger
		return False

	def kGrid(self, startGrid, steps, stepsize = [2,2,2]):
		self.qeio.calculation = 'scf'
		calcType = 'kgrid'
		if self.convEV == 'auto':
			self.convEV == self.qeio.nat*0.001
		self.qeio.kConverged = False
		subDir = self.qeio.kGridFileDir
		opDir = os.path.join(self.qeio.workingDir, subDir)
		if self.replot:
			returnDict = self.qeio.fitKgrid(subDir)
			message = '%s calculation complete on project: %s' % ('k-grid', self.qeio.project)
			self.qeio.report(message, text=self.text, email=self.email, files=returnDict['files'])
			return
		if self.restart: 
			if not os.path.isdir(opDir):
				self.qeio.log('You are attempting to restart a calculation, but the previous results were not found. Set restart=False to start a new calculation.')
				return
			self.qeio.restart = True
		else:
			self.qeio.clearDir(subDir)
		trigger = []
		triggered = False
		converged = False
		self.qeio.kxyz = startGrid
		self.qeio.startingwfc = self.startingwfc
		self.qeio.startingpot = self.startingpot
		if self.convEV == 'auto':
			self.convEV = self.qeio.nat * 0.001
		self.qeio.log('k-grid calculation commenced %s (server time)' % time.asctime())
		for i in range(steps):
			self.qeio.optimization = 'kgrid'
			self.qeio.optLocals = locals()
			opFilename = os.path.join(opDir, '%s-kgrid-%i-%i-%i.out' % (self.qeio.postfix, self.qeio.kxyz[0], self.qeio.kxyz[1], self.qeio.kxyz[2]))#opDir+'/kgrid-%i-%i-%i' % (self.qeio.kxyz[0], self.qeio.kxyz[1], self.qeio.kxyz[2]) 
			self.qeio.write(calcType)
			self.qeio.runPw(opFilename, calcType)
			#self.qeio.startingwfc = 'file'
			#self.qeio.startingpot = 'file'
			data = self.qeio.parseResults(subDir)
			data = sorted(data, key=lambda k: float(float(k['kGridDict']['x']) + float(k['kGridDict']['y']) + float(k['kGridDict']['z'])))
			self.qeio.log('Total energy for %s k-points: %f' % (str(data[-1]['kGridDict']), data[-1]['energy']))
			self.qeio.log('Total time: %is' % data[-1]['wallTime'])
			self.qeio.log('Date/time: %s' % data[-1]['dateTime'])
			if i > 0:
				deltaEV = data[-1]['energy'] - data[-2]['energy']
				self.qeio.log(self.qeio.highlight(str(deltaEV)+' eV change'))
				if abs(deltaEV) < self.convEV:
					if triggered:
						self.qeio.kxyz = trigger
						self.qeio.log('k divisions converged at '+str(self.qeio.kxyz))
						self.qeio.logValue({'k-grid': '%s-%s-%s' % (self.qeio.kxyz[0], self.qeio.kxyz[1], self.qeio.kxyz[2])})
						converged = True
						break
					else:
						triggered = True
						trigger = copy.copy(self.qeio.kxyz)
			self.qeio.kxyz = [self.qeio.kxyz[0] + stepsize[0], self.qeio.kxyz[1] + stepsize[1], self.qeio.kxyz[2] + stepsize[2]]
		self.qeio.parseResults(subDir, log='kGrid')
		returnDict = self.qeio.fitKgrid(subDir, plot=self.plot)
		if self.email or self.text:
			message = '%s calculation complete on project: %s' % ('k-grid', self.qeio.project)
			self.qeio.report(message, text=self.text, email=self.email, files=returnDict['files'])
		if not converged:
			self.qeio.log('k grid did not converge')
			return False
		return trigger

		
	def relax(self, nstep, constraints = False, vc=False, cellDoFree=False):
		self.qeio.calculation = 'scf'
		self.qeio.nstep = nstep
		calcType = 'vc-relax' if vc else 'relax'
		calcTypeString = 'Variable-cell' if vc else 'Fixed-cell'
		subDir = self.qeio.vcRelaxFileDir if vc else self.qeio.relaxFileDir
		self.qeio.constraints = constraints
		self.qeio.calculation = calcType
		if calcType == 'vc-relax' and cellDoFree:
			self.qeio.cellDoFree = cellDoFree
		self.qeio.startingwfc = self.startingwfc
		self.qeio.startingpot = self.startingpot
		self.qeio.write(calcType)
		opDir = os.path.join(self.qeio.workingDir,subDir)
		opFilename = os.path.join(opDir, '%s-%s.out' % (self.qeio.postfix, calcType))
		if not self.replot:
			self.qeio.optimization = 'relax'
			self.qeio.optLocals = locals()
			self.qeio.log('%s relaxation calculation commenced at %s (server time)' % (calcTypeString, time.asctime()))
			if self.restart: 
				if not os.path.isdir(opDir):
					self.qeio.log('You are attempting to restart a calculation, but the previous results were not found. Set restart=False to start a new calculation.')
					return
				self.qeio.restart = True
			else:
				self.qeio.clearDir(subDir)
			self.qeio.runPw(opFilename, calcType)
			self.qeio.parseResults(subDir, log=calcType)
		self.qeio.exportStructure(opFilename, 'espresso-out', self.qeio.postfix+'-'+calcType+'.cif', exportFormat='cif')
		self.qeio.constraints = False
		self.qeio.calculation = 'scf'
		if self.email or self.text:
			message = '%s calculation complete on project: %s' % ('relaxation', self.qeio.project)
			self.qeio.report(message, text=self.text, email=self.email)
		
	def quickRun(self, calcType = 'quickrun', logvbm = False):
		self.qeio.calculation = 'scf'
		subDir = self.qeio.quickRunFileDir
		self.qeio.clearDir(subDir, False)
		fileTail = '%s-%s.out' % (self.qeio.postfix, calcType)
		opFilename = os.path.join(self.qeio.workingDir,subDir,fileTail)
		self.qeio.startingpot = self.startingpot
		self.qeio.startingwfc = self.startingwfc
		self.qeio.write(calcType)
		self.qeio.log('Calculation (%s) commenced %s (server time)' % (calcType, time.asctime()))
		self.qeio.runPw(opFilename, calcType)
		data = self.qeio.parseResults(subDir, files=[fileTail], log=calcType)[0]
		self.qeio.log('Wall time: %is' % data['wallTime'])
		message = '%s calculation complete on project: %s' % (calcType, self.qeio.project)
		self.qeio.report(message, text=self.text, email=self.email)
		if logvbm:
			if data['vbm']:
				self.qeio.vbm = data['vbm']
			if data['cbm']:
				self.qeio.cbm = data['cbm']
			if data['fermiLevel']:
				self.qeio.fermiLevel = data['fermiLevel']
			if data['bandGap']:
				self.qeio.bandGap = data['bandGap']
				self.qeio.logValue({'band-gap':self.qeio.bandGap})

	def bandPlot(self, path, pts, bandLabels = False, nbnd = 'auto', shiftK = False, kxyzMultiplier=2):
		calcType = 'bandplot'
		ogkxyz, ogshiftK = (self.qeio.kxyz, self.qeio.kxyz)
		self.qeio.shiftK = shiftK
		self.qeio.kxyz = np.array(self.qeio.kxyz)*kxyzMultiplier
		self.qeio.wfCollect = True
		if nbnd == 'auto':
			nbnd = self.qeio.valenceElectrons
		self.qeio.nbnd = nbnd
		self.qeio.calculation = 'scf'
		self.quickRun(calcType = 'bandsScf', logvbm=True)
		self.qeio.kxyz = ogkxyz
		self.qeio.kPath = path
		self.qeio.kPathDensity = pts
		subDir = self.qeio.bandPwFileDir
		self.qeio.clearDir(subDir)
		opFilename = os.path.join(self.qeio.workingDir,subDir,'%s-bandplot.out' % self.qeio.postfix)
		self.qeio.calculation = 'bands'
		self.qeio.write(calcType)
		self.qeio.log('Bands non-scf PW calculation commenced at %s (server time)' % time.asctime())
		self.qeio.runPw(opFilename, calcType)
		self.qeio.log('Running bands.x at %s (server time)' % time.asctime())
		self.qeio.runBands(plot = self.plot, bandLabels = bandLabels)
		self.qeio.kPath = False
		self.qeio.shiftK = ogshiftK
		self.qeio.calculation = 'scf'

	def dosPlot(self, Ecenter=False, Espread=False, Emin=False, Emax=False, deltaE = 0.01, nbnd = 'auto', degauss=0.01, kxyzMultiplier=2, noSym=False): # runDOS(self, Ecenter=False, Espread=False, Emin=False, Emax=False, deltaE = 0.01, fildos=False, plot=False):
		calcType, ogkxyz, ogoccupations, ogdegauss, ognoSym = ('dosplot', self.qeio.kxyz, self.qeio.occupations, self.qeio.degauss, self.qeio.noSym)
		self.qeio.occupations = 'smearing'
		self.qeio.degauss = degauss
		self.qeio.wfCollect = True
		if nbnd == 'auto':
			nbnd = self.qeio.valenceElectrons
		self.qeio.nbnd = nbnd
		self.qeio.calculation = 'scf'
		self.quickRun(calcType = 'dosScf', logvbm=True)
		subDir = self.qeio.dosPwFileDir
		self.qeio.clearDir(subDir)
		opFilename = os.path.join(self.qeio.workingDir, subDir, '%s-dosplot.out' % self.qeio.postfix)
		self.qeio.kxyz = np.array(self.qeio.kxyz)*kxyzMultiplier
		self.qeio.calculation = 'nscf'
		self.qeio.write(calcType)
		self.qeio.log('DOS non-scf PW calculation commenced at %s (server time)' % time.asctime())
		self.qeio.runPw(opFilename, calcType)
		result = self.qeio.runDOS(Ecenter=Ecenter, Espread=Espread, Emin=Emin, Emax=Emax, deltaE = deltaE, plot=self.plot)
		self.qeio.occupations = False
		self.qeio.degauss = False
		self.qeio.calculation = 'scf'
		self.qeio.kxyz = ogkxyz
		self.qeio.occupations = ogoccupations
		self.qeio.degauss = ogdegauss
		self.noSym = ognoSym
		if result.errorEncountered:
			if result.errorEncountered == 'dMatrix':
				if self.salvageNumber < 2:
					print('Unable to plot DOS due to d_matrix error')
				else:
					print('Reattempting dosPlot routine due to d_matrix error. Attempting without symmetry.')
					self.salvageNumber += 1
					self.dosPlot(Ecenter=Ecenter, Espread=Espread, Emin=Emin, Emax=Emax, deltaE=deltaE, nbnd=nbnd, degauss=degauss, kxyzMultiplier=kxyzMultiplier, noSym=True)
		else:
			self.salvageNumber = 0


	

		