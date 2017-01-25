import os
import sys
import time
import subprocess
import json
import re
import config
from scipy.optimize import curve_fit
import multiprocessing
import math
if not config.system.hasDisplay:
	import matplotlib
	matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pickle
import optimize
import send
import csv

class io:
	""" This class will write and overwrite an input file for quantum espresso pw.x scf calculations"""
	def __init__(self, project, eCutWfc = 30, eCutRho = 120, a = 6.48, kxyz = [6,6,6], ibrav=2, mpi = False, nbnd=False, cellParams=False):
		self.project = project
		self.bohr2ang = 0.529177
		self.ry2ev = 13.605693
		self.qeDir = config.system.qeDir
		self.pgmDir = config.system.pgmDir
		self.pseudoDir = config.system.pseudoDir
		self.workingDir = os.path.join(self.qeDir, 'projects', self.project)
		self.outdir = os.path.join(self.workingDir, 'results')
		self.logPath = os.path.join(self.workingDir, 'message.log')
		self.prefix = '%s-res' % self.project
		self.mpi = mpi
		self.ibrav = ibrav
		self.cellParams = cellParams
		self.eCutWfc = eCutWfc
		self.eCutRho = eCutRho
		self.elements = {}
		self.atoms = []
		self.ipFile = os.path.join(self.workingDir, 'scf-ip')
		self.kxyz = kxyz
		self.a = a / self.bohr2ang
		self.data = []
		self.calculation = 'scf'
		self.optimization = False
		self.optLocals = False
		self.nstep = False
		self.nat = 0
		self.ntyp = 0
		self.nbnd = nbnd
		self.constraints = False
		self.startingpot = False
		self.startingwfc = False
		self.restart = False
		if not os.path.isdir(self.workingDir):
			try:
				os.mkdir(self.workingDir)
			except:
				print 'Error creating project directory\n'
				print 'Error data:'+sys.exc_info()
		self.rmExitFile()
		if self.ibrav == 0 and not self.cellParams:
			raise Exception("You must manually enter cell tranlation vectors (cellParams) when using the free lattice (ibrav = 0). Set cellParams to a 3 element list of three element cartesian translation vectors (floats). ")
	def addElement(self,symbol,mass,pseudoFile):
		"""Add a reference for an element and it's corresponsing pseudopotential file. Mass in A.U. must also be entered"""
		self.elements[symbol] = {}
		self.elements[symbol]['mass'] = mass
		self.elements[symbol]['pseudoFile'] = pseudoFile
		if not os.path.isfile(os.path.join(self.pseudoDir,pseudoFile)):
			print 'Pseudopotential file not located. There must be a pseudopotential file of UPF format located in the Quantum Espresso pseudo directory'
			return False
		self.renumber()
	def addAtom(self,symbol,position):
		"""Adds a single atom of defined by symbol(element) and x,y,z position list in units of self.a"""
		if symbol not in self.elements:
			print 'That element has not yet been defined. Define elements using the addElement method'
			return False
		self.atoms.append([symbol,position[0],position[1],position[2]])
		self.renumber()
	def importAtomFile(self,filename):
		"""Clears all atoms and add atoms according to the filename file with a json format corresponding to an array of arrays, each with 4 elements like [element symbol, x, y, z] where x,y,z are in units of self.a"""
		if not os.path.isfile(filename):
			print 'Cannot find atoms file.'
			return False
		self.atoms = json.loads(open(filename,'r').read())
		self.renumber()
	def exportAtoms(self, dir = 'atomfiles', filename = 'atoms.atoms'):
		if len(self.atoms) > 0:
			opDir = os.path.join(self.workingDir,dir)
			filepath = os.path.join(opDir,filename)
			if not os.path.isdir(opDir):
				try:
					os.mkdir(opDir)
				except:
					print 'Error creating project directory\n'
					print 'Error data:'+str(sys.exc_info())
					return
			file = open(filepath,'w+')
			file.truncate()
			file.write(json.dumps(self.atoms))
			file.close()
	def renumber(self):
			typArr = []
			tmpNat = 0
			tmpNtyp = 0
			for atom in self.atoms:
				tmpNat += 1
				if atom[0] not in typArr:
					typArr.append(atom[0])
					tmpNtyp += 1
			self.nat = tmpNat
			self.ntyp = tmpNtyp
	def write(self):
		inputFile = open(self.ipFile,'w+')
		inputFile.truncate()
		inputFile.write(" &CONTROL\n")
		inputFile.write("  calculation = '%s',\n" % self.calculation)
		if self.restart:
			inputFile.write("  restart_mode = 'restart',\n")
		inputFile.write("  outdir = '%s',\n" % self.outdir)
		inputFile.write("  pseudo_dir = '%s',\n" % self.pseudoDir)
		inputFile.write("  prefix = '%s',\n" % self.prefix)
		if self.nstep:
			inputFile.write("  nstep = %i,\n" % self.nstep)
		inputFile.write(" /\n")
		inputFile.write(" &SYSTEM\n")
		inputFile.write("  ibrav = %i,\n" % self.ibrav)
		inputFile.write("  celldm(1) = %f,\n" % self.a)
		inputFile.write("  nat = %i,\n" % self.nat)
		inputFile.write("  ntyp = %i,\n" % self.ntyp)
		inputFile.write("  ecutwfc = %d,\n" % self.eCutWfc)
		inputFile.write("  ecutrho = %d\n" % self.eCutRho)
		inputFile.write(" /\n")
		inputFile.write(" &ELECTRONS\n")
		if self.startingwfc:
			inputFile.write("  startingwfc = '%s'\n" % self.startingwfc)
		if self.startingpot:
			inputFile.write("  startingpot = '%s'\n" % self.startingpot)
		inputFile.write(" /\n")
		if self.calculation == 'relax' or 'vc-relax':
			inputFile.write(" &IONS\n")
			if self.constraints and self.calculation == 'relax':
				inputFile.write(" ion_dynamics = 'damp'\n")
			inputFile.write(" /\n")
			if self.calculation == 'vc-relax':
				inputFile.write(" &CELL\n")
				inputFile.write(" /\n")
		inputFile.write("ATOMIC_SPECIES\n")
		for symbol,data in self.elements.iteritems():
			inputFile.write(" %s %f %s\n" % (symbol, data['mass'], data['pseudoFile']))
		inputFile.write("ATOMIC_POSITIONS alat\n")
		for atom in self.atoms:
			inputFile.write(" %s %f %f %f\n" % (atom[0], atom[1], atom[2], atom[3]))
		inputFile.write("K_POINTS automatic\n")
		inputFile.write("%d %d %d 0 0 0\n" % (self.kxyz[0], self.kxyz[1], self.kxyz[2]))
		if self.ibrav == 0:
			inputFile.write("CELL_PARAMETERS\n")
			inputFile.write(" %f %f %f\n" % (self.cellParams[0][0],  self.cellParams[0][1], self.cellParams[0][2]))
			inputFile.write(" %f %f %f\n" % (self.cellParams[1][0],  self.cellParams[1][1], self.cellParams[1][2]))
			inputFile.write(" %f %f %f\n" % (self.cellParams[2][0],  self.cellParams[2][1], self.cellParams[2][2]))
		if (self.calculation == 'relax' or self.calculation == 'vc-relax') and self.constraints:
			inputFile.write("CONSTRAINTS\n")
			if self.constraints == 'zb-bennett':
				inputFile.write("1\n")
				inputFile.write("'bennett_proj' 0 1.0 1.0 1.0 \n")
		inputFile.close()
	def runPw(self, opFilename):
		"""Runs pw.x with the current input file and writes the output to opFilename. self.mpi is a 5-element array that reads [np,ni,nk,nb,nt,nd] (see QE user guide page 23). if self.mpi is set, the calculation will be run in parallel using mpirun"""
		pickleFile = open(os.path.join(self.workingDir, 'io.pkl'),'w')
		pickleFile.truncate()
		pickle.dump(self, pickleFile, pickle.HIGHEST_PROTOCOL)
		if self.mpi:
			nprocs = multiprocessing.cpu_count()
			if len(self.mpi) < 6:
				print 'Your mpi parameters array is too short'
				exit()
			if self.mpi[0] < nprocs:
				print 'Warning: you are running in parallel, but not utilizing all of your processors'
			if self.mpi[0] > nprocs:
				print 'You have specified more processors than are available'
				exit()
			if self.mpi[1]*self.mpi[2] > self.mpi[0]:
				print 'Your mpi paramters specify too many images/pools. Use constraint ni*nk <= np'
				exit()
			if self.mpi[5] > self.mpi[0]/(self.mpi[1]*self.mpi[2]):
				print 'You have specified too many processors for your diagonalization group.'
				exit()
			if not self.mpi[5] == int(math.sqrt(self.mpi[5]))**2:
				print 'The diagonalization group size must be a perfect square'
				exit()
			if not isinstance(self.mpi[2]/self.mpi[4], int):
				print 'Warning: your mpi params are not optimized. Try constraint nk/nt = Int'
			cmdStr = ' '.join(['mpirun', '-np', str(self.mpi[0]), os.path.join(self.pgmDir,'pw.x'), '-ni', str(self.mpi[1]), '-nk', str(self.mpi[2]), '-nb', str(self.mpi[3]), '-nt', str(self.mpi[4]), '-nd', str(self.mpi[5]), '-i', self.ipFile, '>', opFilename])
			os.system(cmdStr)
			self.appendData(opFilename)
		else:
			os.system("%s < %s > %s" % (self.pgmDir+'/pw.x', self.ipFile, opFilename))
			self.appendData(opFilename)
	def parseResults(self, dir, log=False):
		"""Reads output files from an optimization, and creates an array of data"""
		fileDir = self.sub2path(dir)
		if not os.path.isdir(fileDir):
			print 'Could not parse results. Directory not found. Directory: %s' % fileDir
			return False
		self.data = []
		reLattice = re.compile('lattice parameter')
		reEnergy = re.compile('!')
		reCutWfc = re.compile('kinetic-energy cutoff')
		reCutRho = re.compile('charge density cutoff')
		reTime = re.compile(r'PWSCF[\s]*:[ \w.]* WALL')
		reNonDec = re.compile(r'[^\d.]+')
		reKpts = re.compile('number of k points')
		reKgrid = re.compile('kGrid')
		reMpi = re.compile('mpi:')
		reProj = re.compile('project:')
		reDate = re.compile('Local time:')
		reEpoch = re.compile('epochTime:')
		reStartingwfc = re.compile('startingwfc:')
		reStartingpot = re.compile('startingpot:')
		for filename in os.listdir(fileDir):
			pt = {}
			file = open(os.path.join(fileDir,filename))
			for line in file.readlines():
				if(re.search(reLattice, line)):
					words = line.split()
					pt['latticeConstant'] = float(words[-2])*self.bohr2ang
				if(re.search(reEnergy, line)):
					words = line.split()
					pt['energy'] = float(words[-2])*self.ry2ev
				if(re.search(reCutWfc, line)):
					words = line.split()
					pt['ecutwfc'] = float(words[-2])
				if(re.search(reTime, line)):
					try:
						split1 = line.split(':')
						split2 = split1[1].strip().split('CPU')
						pt['cpuTime'] = self.str2secs(split2[0].strip())
						split3 = split2[1].strip().split('WALL')
						pt['wallTime'] = self.str2secs(split3[0].strip())
					except:
						print 'Error parsing time from line: '+line
						print 'Error data:'+str(sys.exc_info())
						pt['wallTime'] = 1
						pt['cpuTime'] = 1
				if(re.search(reCutRho, line)):
					words = line.split()
					pt['ecutrho'] = float(words[-2])
				if(re.search(reKpts, line)):
					words = line.split()
					pt['kPts'] = float(words[-1])
				if(re.search(reKgrid, line)):
					matches = re.search('\[(?P<x>\d+), (?P<y>\d+), (?P<z>\d+)\]', line)
					pt['kGridDict'] = matches.groupdict()
				if(re.search(reMpi, line)):
					matches = re.search('\[(?P<np>\d+), (?P<ni>\d+), (?P<nk>\d+), (?P<nb>\d+), (?P<nt>\d+), (?P<nd>\d+)\]', line)
					pt['mpi'] = matches.groupdict()
				if(re.search(reProj, line)):
					pt['project'] = line.strip().split(':')[1]
				if(re.search(reDate, line)):
					pt['dateTime'] = line.strip().split(':', 1)[1]
				if(re.search(reEpoch, line)):
					pt['epochTime'] = float(line.strip().split(':', 1)[1])
				if(re.search(reStartingwfc, line)):
					pt['startingwfc'] = line.strip().split(':', 1)[1]
				if(re.search(reStartingpot, line)):
					pt['startingpot'] = line.strip().split(':', 1)[1]
			if log:
				pt['calcType'] = log
			pt['nat'] = self.nat
			self.data.append(pt)
		if log:
			filePath = os.path.join(self.workingDir, 'data.log')
			if os.path.isfile(filePath):
				with open(filePath, 'r') as file:
					dataList = json.loads(file.read())
			else:
				dataList = []
			dataList = dataList + self.data
			with open(filePath, 'w+') as file:
				file.write(json.dumps(dataList))
		return self.data
	def fitLattice(self, plot=False, subDir='lattice-op'):
		""" Plots the output of self.parseResults after a lattice optimization """
		self.parseResults(subDir)
		opDir = self.sub2path(subDir)
		x = []
		y = []
		pts = sorted(self.data, key=lambda k: float(k['latticeConstant']))
		if len(pts) < 5:
			print "You don't appear to have enough data points. Be sure to collect at least 5 data points and to call qeio.parseResults before fitting the lattice data."
			return False
		for pt in pts:
			y.append(pt['energy'])
			x.append(pt['latticeConstant'])
		def quad(a,E0,B,a0):
			return E0 + B*(a - a0)**2
		quadParams, pcov = curve_fit(quad, x, y, p0=[y[1],1,x[1]])
		print str(quadParams)
		def BMeos(a,B0,a0,dB,E0,test):
			A = (a0/a)**2
			return E0 + test*((A - 1)**3*dB + (A - 1)**2*(6 - 4*A))
		BMparams, pcov = curve_fit(BMeos, x, y, p0=[1,x[1],1,y[1],1])
		if plot:
			denseX = []
			tmpX = min(x)
			deltaX = (max(x) - min(x))/30
			while tmpX < max(x):
				denseX.append(tmpX)
				tmpX += deltaX
			fig = plt.figure()
			fig.patch.set_facecolor('white')			
			dataMarkers = plt.plot(x,y)
			plt.setp(dataMarkers, color='k', marker='s', linestyle = 'none', markersize = 8.0, label = 'Data')
			quadY = [quad(a, quadParams[0], quadParams[1], quadParams[2]) for a in denseX]
			quadLine = plt.plot(denseX, quadY)
			plt.setp(quadLine,color='k', marker=None, linestyle='solid', linewidth=1.0, label = 'Quadratic')
			BMy = [BMeos(a, BMparams[0], BMparams[1], BMparams[2], BMparams[3], BMparams[4]) for a in denseX]
			BMline = plt.plot(denseX, BMy)
			plt.setp(BMline,color='k', marker=None, linestyle='dotted', linewidth=1.0, label = 'BM')
			plt.legend(loc="upper center", fancybox=True, ncol=1, fontsize=13, numpoints = 1)
			plt.suptitle('DFT energy v. Lattice Constant for FCC CdTe', fontsize=20)
			plt.xlabel(r'Lattice constant ($\AA$)', fontsize=18)
			plt.ylabel('DFT energy (Ry)', fontsize=18)
			plt.subplots_adjust(left = 0.2, right = 0.9, top = 0.85, bottom = 0.2)
			filePath = '.'.join([os.path.join(self.workingDir, subDir), 'png'])
			fig.savefig(filePath,format='png')
			returnDict = {}
			returnDict['files'] = []
			returnDict['files'].append(filePath)
			returnDict['BMfit'] = BMparams[1]
			returnDict['quadFit'] = quadParams[2]
			if config.system.hasDisplay:
				plt.ion()
				plt.show()
			else:
				plt.close()
		return returnDict
	def fitCutoff(self, qeParam, subDir):
		""" Plots and fits trendlines to data from energy cutoff convergence calculation"""
		self.parseResults(subDir)
		x = []
		y = []
		t = []
		pts = sorted(self.data, key=lambda k: float(k[qeParam]))
		if len(pts) < 5:
			print "You don't appear to have enough data points. Be sure to collect at least 5 data points and remember to call qeio.parseResults before fitting the data."
			return False
		for pt in pts:
			y.append(pt['energy'])
			print 'energy from point %i: %f' % (len(y), pt['energy'])
			if len(y) > 1:
				print 'energy change: %s' % str(y[-1] - y[-2])
			x.append(pt[qeParam])
			t.append(pt['wallTime'])
		return self.fitConvergence(x, y, t, subDir, xLabel=qeParam+' Cutoff')
	def fitKgrid(self, subDir, plot=True):
		""" Plots and fits trendlines to data from density energy cutoff convergence calculation"""
		self.parseResults(subDir)
		x = []
		y = []
		t = []
		xTickLabels = []
		pts = sorted(self.data, key=lambda k: float(k['kPts']))
		if len(pts) < 5:
			print "You don't appear to have enough data points. Be sure to collect at least 5 data points and to call qeio.parseResults before fitting the data."
			return False
		for pt in pts:
			y.append(pt['energy'])
			x.append(float(pt['kGridDict']['x']) + float(pt['kGridDict']['y']) + float(pt['kGridDict']['z']))
			xTickLabels.append('[%s,%s,%s]' % (pt['kGridDict']['x'], pt['kGridDict']['y'], pt['kGridDict']['z']))
			t.append(pt['wallTime'])
		return self.fitConvergence(x, y, t, subDir, xLabel = 'k-Grid Divisions', xTickLabels=xTickLabels)
	def fitConvergence(self, x, y, t, subDir, xLabel='', xTickLabels=False):
		deltaX = x[1:]
		deltaY = [abs(j-i) for i,j in zip(y[:-1], y[1:])]
		denseX = []
		tmpX = min(x)
		xStep = (max(x) - min(x))/30
		while tmpX < max(x):
			denseX.append(tmpX)
			tmpX += xStep
		fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)
		plt.suptitle('Convergence Plot', fontsize=20)
		plt.xlabel(xLabel, fontsize=18)
		if xTickLabels:
			plt.xticks(x, xTickLabels)
		plt.subplots_adjust(left = 0.2)
		fig.patch.set_facecolor('white')
		fig.xlabel = 'Cutoff Energy (eV)'		
		ax1.plot(x, y, marker='s', color='k', linestyle = 'none', markersize = 5.0, label = 'Data')
		ax1.set_ylabel('DFT energy (eV)')
		ax2.plot(deltaX, deltaY, color='k', marker='s', linestyle = 'none', markersize = 5.0, label = 'Data')
		ax2.set_ylabel(r'$\Delta$E (eV)')
		ax3.plot(x, t, marker='s', color='k', linestyle = 'none', markersize = 5.0, label = 'Data')
		ax3.set_ylabel('Time (s)')
		filePath = '.'.join([os.path.join(self.workingDir, subDir), 'png'])
		fig.savefig(filePath,format='png')
		returnDict = {}
		returnDict['files'] = []
		returnDict['files'].append(filePath)
		if config.system.hasDisplay:
			plt.ion()
			plt.show()
			raw_input('Press Enter button to continue')
		else:
			plt.close()
		return returnDict
	def xcrysden(self, filename = False):
		if not config.system.hasDisplay:
			print 'xcrysden needs a display to run. Your configuration file indicates that you have no display. '
			return
		self.write()
		if not filename:
			filepath = self.ipFile
		if config.system.hasDisplay:
			subprocess.call('xcrysden --pwi '+filepath, shell = True)
	def str2secs(self, string):
		tString = string.strip()
		t=0
		match = re.search('\D*(\d+)d',tString)
		if match:
			t += 86400*float(match.group(1))
		match = re.search('\D*(\d+)h',tString)
		if match:
			t += 3600*float(match.group(1))
		match = re.search('\D*(\d+)m',tString)
		if match:
			t += 60*float(match.group(1))
		match = re.search('\D*([.\d]+)s',tString)
		if match:
			t += float(match.group(1))
		return t
	def clearDir(self, subDir):
		""" Clears directory. Creates empty directory if it doesn't exist"""
		dirPath = os.path.join(self.workingDir, subDir)
		if not os.path.isdir(dirPath):
			try:
				os.mkdir(dirPath)
			except:
				print 'Error creating output directory\n'
				print 'Error data:'+str(sys.exc_info())
				return False
		for file in os.listdir(dirPath):
			filePath = os.path.join(dirPath, file)
			try:
				if os.path.isfile(filePath):
					os.unlink(filePath)
			except:
				print 'Error clearing lattice output directory'
				print 'Error data:'+sys.exc_info()
				return False
	def sub2path(self, subDir):
		return os.path.join(self.workingDir, subDir)
	def appendData(self, filepath):
		with open(filepath, 'a') as file:
			file.write('\n kGrid:%s' % str(self.kxyz))
			file.write('\n project:%s' % self.project)	
			if self.mpi:
				file.write('\n mpi:%s' % str(self.mpi))
			else:
				file.write('\n mpi:[1, 1, 1, 1, 1, 1]')
			file.write('\n Local time: %s' % time.asctime())
			file.write('\n startingwfc: %s' % self.startingwfc)
			file.write('\n startingpot: %s' % self.startingpot)
			file.write('\n epochTime: %s' % time.time())
	def rmExitFile(self):
		filepath = os.path.join(self.outdir,self.prefix+'.EXIT')
		if os.path.isfile(filepath):
			os.remove(filepath)
	def pause(self):
		if os.path.isdir(self.outdir):
			filepath = os.path.join(self.outdir,self.prefix+'.EXIT')
			open(filepath, 'a').close()
			print 'Kill switch activated. If any calculations are running for the provided project, they should end shortly.'
		else:
			print 'There does not seem to be any project with that name.'
	def restartCalc(self):
		if self.optimization == 'lattice':
			optimize.latticeConstant(self, low = self.optLocals['low'], high = self.optLocals['high'], steps = self.optLocals['steps'], plot = self.optLocals['plot'], replot = self.optLocals['replot'], restart = True)
		elif self.optimization == 'ecutwfc':
			optimize.cutWfc(self, low = self.optLocals['testE'], high = self.optLocals['high'], stepsize = self.optLocals['stepsize'], convEV = self.optLocals['convEV'], plot = self.optLocals['plot'], replot = self.optLocals['replot'], restart=True)
		elif self.optimization == 'ecutrho':
			optimize.cutRho(self, low = self.optLocals['testE'], high = self.optLocals['high'], stepsize = self.optLocals['stepsize'], convEV = self.optLocals['convEV'], plot = self.optLocals['plot'], replot = self.optLocals['replot'], restart=True)
		elif self.optimization == 'kgrid':
			optimize.kGrid(self, low = self.optLocals['testDivs'], high = self.optLocals['high'], convEV = self.optLocals['convEV'], plot = self.optLocals['plot'], replot = self.optLocals['replot'], restart=True)
		elif self.optimization == 'relax':
			optimize.relax(self, nstep = self.optLocals['nstep'], constraints = self.optLocals['constraints'], vc = self.optLocals['vc'], replot = self.optLocals['replot'], restart=True)
		else:
			print "There doesn't appear to be a calculation to restart"
	def log(self, entry):
		if not entry:
			return
		print entry
		with open(self.logPath, 'a') as logFile:
			logFile.write(str(entry)+'\n')
	def report(self, message, email=False, subject = '', text=False, files=False, image=False):
		if email:
			send.email(subject, message, files=files)
		if text:
			send.text(message)
	def exportData(self):
		importPath = os.path.join(self.workingDir, 'data.log')
		exportPath = os.path.join(self.workingDir,self.project+'.csv')
		if os.path.isfile(importPath):
			with open(importPath, 'r') as file:
				dataList = json.loads(file.read())
			with open(exportPath, 'w+') as file:
				writer = csv.DictWriter(file, self.flattenData(dataList[0]).keys())
				writer.writeheader()
				for dataPt in dataList:
					writer.writerow(self.flattenData(dataPt))
		else:
			print 'No data to export.'
	def flattenData(self, dataPt):
		opDict = {}
		if 'mpi' not in dataPt:
			dataPt['mpi'] = {'np':'1','ni':'1','nk':'1','nb':'1','nt':'1','nd':'1'}
		for key in dataPt:
			if hasattr(dataPt[key], '__iter__'):
				for subkey in dataPt[key]:
					opDict['.'.join([key,subkey])] = dataPt[key][subkey]
			else:
				opDict[key] = dataPt[key]			
		return opDict
			
			