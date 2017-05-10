from compatibility import *
import os
import shutil
import sys
import time
from subprocess import Popen, PIPE
import json
import re
import config
from scipy.optimize import curve_fit
import numpy.linalg as la
import numpy as np
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
import pointgroup as pg
from collections import OrderedDict

class GenericClass(object):
	pass

class io:
	""" This class will write and overwrite an input file for quantum espresso pw.x scf calculations"""
	def __init__(self, project, pseudoType = False, ecutwfc = 30, ecutrho = 120, a = 6.48, kxyz = [6,6,6], ibrav=0, mpi = False, cellParams=[]):
		self.project = project
		self.bohr2ang = 0.529177249
		self.ry2ev = 13.605693
		self.qeDir = config.system.qeDir
		self.pgmDir = config.system.pgmDir
		self.pseudoDir = config.system.pseudoDir
		self.structureDir = config.system.structureDir
		self.projectDir = config.system.projectDir
		self.dosFileDir = 'dos-files'
		self.bandFileDir = 'band-files'
		self.latticeFileDir = 'lattice-files'
		self.ecutwfcFileDir = 'ecutwfc-files'
		self.ecutrhoFileDir = 'ecutrho-files'
		self.kGridFileDir = 'kgrid-files'
		self.vcRelaxFileDir = 'vcrelax-files'
		self.relaxFileDir = 'relax-files'
		self.quickRunFileDir = 'quickrun-files'
		self.bandPwFileDir = 'band-pw-output'
		self.dosPwFileDir = 'dos-pw-output'
		self.inputFileFolder = 'input-files'
		self.defaultNames()
		self.pseudoChecked = False
		self.pseudoType = pseudoType
		self.pseudoFiles = {}
		self.relaxedStructure = False # Hold the filename of the last exported structure
		self.exportFormat = 'cif' # Format of last exported structure file
		self.mpi = mpi
		self.lastError = False
		self.waitAtImages = True
		self.ibrav = ibrav
		self.verbosity = 'low'
		self.wfCollect = False
		self.ecutwfc = ecutwfc
		self.ecutrho = ecutrho
		self.atoms = []
		self.crystal = pg.Crystal(self)
		for v in cellParams:
			self.crystal.addTranslation(v)
		self.kxyz = self.ogkxyz = kxyz
		self.shiftK = False
		self.kPath = False
		self.kPathDensity = False
		self.crystal.a = a 
		self.data = []
		self.calculation = 'scf'
		self.cellDoFree = False
		self.noSym = False
		self.optimization = False
		self.optLocals = False
		self.nstep = False
		self.nat = 0
		self.ntyp = 0
		self.valenceElectrons = False # will be set every time self.parseResults is run. Can be used for band and DOS plotting in nbnd='auto' mode
		self.nbnd = False
		self.bandMultiplier = 1
		self.occupations = False
		self.degauss = False
		self.input_dft = False
		self.fockMesh = False
		self.hubbardU = False # syntax: dictionary with key as atom types and values as U values in eV
		self.ldaUtype1 = False
		self.charge = False
		self.relativistic = False
		self.lspinorb = False
		self.constraints = False
		self.startingpot = False
		self.startingwfc = False
		self.restart = False
		self.dosEmin = False
		self.dosEmax = False
		self.dosDeltaE = 0.01
		self.dosEcenter = False
		self.vbm = False
		self.cbm = False
		self.fermiLevel = False
		self.bandGap = False
		self.rmExitFile()
		self.initDirs()
		self.optimize = optimize.Optimizer(self)
		self.salvageNumber = 0
	def reset(self):
		"""Set many optional values back to default values"""
		self.kxyz = self.ogkxyz
		self.shiftK = False
		self.kPath = False
		self.kPathDensity = False
		self.cellDoFree = False
		self.noSym = False
		self.nbnd = False
		self.bandMultiplier = 1
		self.occupations = False
		self.degauss = False
		self.input_dft = False
		self.fockMesh = False
		self.hubbardU = False
		self.charge = False
		self.constraints = False
		self.vbm = False
		self.cbm = False
		self.fermiLevel = False
		self.bandGap = False
		self.relativistic = False
		self.lspinorb = False
	def defaultNames(self):
		"""Sets working directories and filenames to default values. May need to call after subtests to reset directory names"""
		self.workingDir = os.path.join(self.projectDir, self.project)
		self.prefix = '%s-res' % self.project
		self.postfix = self.project		
		self.outdir = os.path.join(self.workingDir, 'results')
		self.logPath = os.path.join(self.workingDir, '%s-message.log' % self.postfix)
		self.dictionaryPath = os.path.join(self.workingDir, '%s.dictionary' % self.postfix)
		self.ipFileDir = os.path.join(self.workingDir, self.inputFileFolder)
		self.ipDOS = os.path.join(self.ipFileDir, '%s-dosx.in' % self.postfix)
		self.ipBands = os.path.join(self.ipFileDir, '%s-bandsx.in' % self.postfix)
		self.dosPath = os.path.join(self.workingDir,self.dosFileDir)
		self.bandPath = os.path.join(self.workingDir, self.bandFileDir)
		self.fildos = os.path.join(self.dosPath, self.postfix)
		self.filband = os.path.join(self.bandPath, '%s.bands' % self.postfix)
		self.datalog = os.path.join(self.workingDir, '%s-data.log' % self.postfix)
		self.csvPath = os.path.join(self.workingDir, '%s-data.csv' % self.postfix)
	def subTest(self, subtest = False, chain = False):
		"""Creates and works in a subdirectory in the working directory (self.workingDir). If called without an argument, or with an argument of False, it will revert directory names to defaults"""
		if subtest:
			if not chain:
				self.defaultNames()
			self.workingDir = os.path.join(self.workingDir, subtest)
			self.prefix = '%s-%s' % (subtest,self.prefix)
			self.postfix = '%s-%s' % (self.postfix, subtest)
			self.outdir = os.path.join(self.workingDir, 'results')
			self.logPath = os.path.join(self.workingDir, '%s-message.log' % self.postfix)
			self.dictionaryPath = os.path.join(self.workingDir, '%s.dictionary' % self.postfix)
			self.ipFileDir = os.path.join(self.workingDir, self.inputFileFolder)
			self.ipDOS = os.path.join(self.ipFileDir, '%s-dosx.in' % self.postfix)
			self.ipBands = os.path.join(self.ipFileDir, '%s-bandsx.in' % self.postfix)
			self.dosPath = os.path.join(self.workingDir,self.dosFileDir)
			self.bandPath = os.path.join(self.workingDir, self.bandFileDir)
			self.fildos = os.path.join(self.dosPath, self.postfix)
			self.filband = os.path.join(self.bandPath, '%s.bands' % self.postfix)
			self.datalog = os.path.join(self.workingDir, '%s-data.log' % self.postfix)
			self.csvPath = os.path.join(self.workingDir, '%s-data.csv' % self.postfix)
			self.initDirs()
		else:
			self.defaultNames()
	'''
	def addElement(self,symbol,mass=False,pseudoFile=False):
		"""Add a reference for an element and it's corresponsing pseudopotential file. Mass in A.U. must also be entered"""
		self.crystal.addElement(symbol,mass,pseudoFile)
		if pseudoFile:
			if not os.path.isfile(os.path.join(self.pseudoDir,pseudoFile)):
					print('Pseudopotential file not located. There must be a pseudopotential file of UPF format located in the Quantum Espresso pseudo directory')
					return False
		self.renumber()
	'''
	def setPseudo(self, symbol, filepath):
		"""Create an entry in the self.pseudoFiles for the given element"""
		if os.path.isfile(os.path.join(config.system.pseudoDir, filepath)):
			self.pseudoFiles[symbol] = filepath
			return True
		return False
	def setPseudoType(self, type):
		""" Set the pseudopotential type. Must match a pseudopotential file in the config.system.pseudoDir directory with the format [atomic symbol].[type].UPF"""
		self.pseudoType = type
		self.pseudoChecked = False
		self.pseudoFiles = {}
		self.checkPseudos()
	def checkPseudos(self):
		"""Check for the existence of pseudopotential files for all elements in the crystal structure"""
		for symbol in self.crystal.elements:
			if(not self.checkPseudo(symbol)):
				exit('Pseudopotential file of type %s for element %s not found. ' % (self.pseudoType, symbol))
		self.pseudoChecked = True
		return True
	def checkPseudo(self, symbol, filename = False):
		"""Check for the existence of pseudopotential file for given element pseudopotential type"""
		if not filename:
			filename = '%s.%s.(upf|UPF)' % (symbol, self.pseudoType)
		pseudoList = os.listdir(config.system.pseudoDir)
		listString = ' '.join(pseudoList)
		m = re.search(filename, listString)
		if m:
			self.pseudoFiles[symbol] = m.group(0)
			return True
		return False
	def fetchPseudo(self, symbol):
		""" Return filename defined by format [symbol].[self.pseudoType].UPF is in directory defined at config.system.pseudoDir"""
		filename = '%s.%s.UPF' % (symbol, self.pseudoType)
		return filename if os.path.isfile(os.path.join(config.system.pseudoDir, filename)) else False
	def makeRelativistic(self):
		""" For fully relativistic pseudopotentials, adjust input file parameters as necessary"""
		self.relativistic = True
		self.lspinorb = True
		self.bandMultiplier = 2
		self.ldaUtype1 = True
	def setA(self, a):
		self.crystal.a = a
	def addAtom(self,symbol,position,forceX=[1,1,1]):
		"""Adds a single atom of defined by symbol(element) and x,y,z position list in units of self.crystal.a"""
		self.crystal.addAtom(symbol, position, forceX)
		self.renumber()
	def removeAtom(self, element, position):
		"""Remove atom from self.crystal. Passes to pointgroup.Crystal::removeAtom. Here as a convenience. """
		return self.crystal.removeAtom(element, position)
	def replaceAtom(self, element, replacement):
		"""Passes to pointgroup.Crystal::replaceAtom. Here as a convenience. Remove the first atom in self.crystal.basis of type element and replace with an atom of type replacement"""
		self.crystal.replaceAtom(element, replacement)
		self.renumber()
	def addAtoms(self, atoms):
		for atom in atoms:
			self.addAtom(atom[0],atom[1],atom[2])
	def repeatStructure(self, t, n):
		"""Passes to pointgroup.Crystal::repeatStructure. Here as a convenience. Repeat the given structure n times along translation vector t, and adjust translation vector t as necessary for supercell."""
		self.crystal.repeatStructure(t, n);
	def printStructure(self):
		"""Print a list of atoms in self.crystal.basis. Passes to PointGroup.Crystal::printStructure. Here as a convenience"""
		self.crystal.printStructure()
	def importStructure(self, path, directory=False, format=None, index=None):
		"""Clears all atoms and add atoms according to the filename file with a json format corresponding to an array of arrays, each with 4 elements like [element symbol, x, y, z] where x,y,z are in units of self.crystal.a"""
		self.crystal.cleanStart()
		directory = directory if directory else self.structureDir
		path = os.path.join(directory,path)
		if not os.path.isfile(path):
			print('Cannot find structure file.')
			return False
		if format == 'json':
			structure = json.loads(open(path,'r').read())
			self.crystal.basis = structure['basis']
			self.crystal.translations = structure['translations']
		else:
			self.crystal.importStructure(path,format=format, index=index)
			self.setA(self.crystal.a)
		self.renumber()
	def exportStructure(self, importPath, importFormat, exportName, exportFormat='cif'):
		"""Export cell and atom info. importFormat will be either espresso-out or espresso-in"""
		self.relaxedStructure = exportName
		self.exportFormat = exportFormat
		exportPath = os.path.join(self.structureDir,exportName)
		self.crystal.exportStructure(importPath, importFormat=importFormat, index=None, exportPath=exportPath, exportFormat=exportFormat)
	def loadRelaxedStructure(self):
		self.importStructure(self.relaxedStructure, format = self.exportFormat)
	def renumber(self):
			typArr = []
			tmpNat = 0
			tmpNtyp = 0
			for atom in self.crystal.basis:
				tmpNat += 1
				if atom[0] not in typArr:
					typArr.append(atom[0])
					tmpNtyp += 1
			self.nat = tmpNat
			self.ntyp = tmpNtyp
	def parseResults(self, subDir, log=False, files=False):
		"""Reads output files from an optimization, and creates an array of data"""
		fileDir = self.sub2path(subDir)
		if not os.path.isdir(fileDir):
			print('Could not parse results. Directory not found. Directory: %s' % fileDir)
			return False
		self.data = []
		if files:
			filePaths = [os.path.join(fileDir, i) for i in files]
		else:
			filePaths = os.listdir(fileDir)
		for filename in filePaths:
			if not filename.endswith('.out'):
				continue
			pt = {}
			pt['latticeConstant'] = ''
			pt['energy'] = ''
			pt['ecutwfc'] = ''
			pt['cpuTime'] = ''
			pt['wallTime'] = ''
			pt['ecutrho'] = ''
			pt['kPts'] = ''
			pt['kGridDict'] = ''
			pt['mpi'] = ''
			pt['project'] = ''
			pt['dateTime'] = ''
			pt['epochTime'] = ''
			pt['startingwfc'] = ''
			pt['startingpot'] = ''
			pt['fermiLevel'] = ''
			pt['cbm'] = ''
			pt['vbm'] = ''
			pt['bandGap'] = ''
			pt['calcType'] = ''
			pt['V'] = ''
			pt['XCfunc'] = ''
			pt['valenceElectrons'] = ''
			file = open(os.path.join(fileDir,filename))
			for line in file.readlines():
				if(re.search('lattice parameter', line)):
					words = line.split()
					pt['latticeConstant'] = float(words[-2])*self.bohr2ang
				if(re.search('!', line)):
					words = line.split()
					pt['energy'] = float(words[-2])*self.ry2ev
				if(re.search('kinetic-energy cutoff', line)):
					words = line.split()
					pt['ecutwfc'] = float(words[-2])
				if(re.search(r'PWSCF[\s]*:[ \w.]* WALL', line)):
					try:
						split1 = line.split(':')
						split2 = split1[1].strip().split('CPU')
						pt['cpuTime'] = self.str2secs(split2[0].strip())
						split3 = split2[1].strip().split('WALL')
						pt['wallTime'] = self.str2secs(split3[0].strip())
					except:
						print('Error parsing time from line: '+line)
						print('Error data:'+str(sys.exc_info()))
						pt['wallTime'] = 1
						pt['cpuTime'] = 1
				if(re.search('charge density cutoff', line)):
					words = line.split()
					pt['ecutrho'] = float(words[-2])
				if(re.search('number of k points', line)):
					m = re.search(r"number of k points=[ ]*([0-9.]+)[\D ]*", line)
					pt['kPts'] = float(m.group(1))
				if(re.search('kGrid', line)):
					matches = re.search('\[(?P<x>\d+)[, ]*(?P<y>\d+)[, ]*(?P<z>\d+)\]', line)
					pt['kGridDict'] = matches.groupdict()
				if(re.search('mpi:', line)):
					matches = re.search('\[(?P<np>\d+), (?P<ni>\d+), (?P<nk>\d+), (?P<nb>\d+), (?P<nt>\d+), (?P<nd>\d+)\]', line)
					pt['mpi'] = matches.groupdict()
				if(re.search('project:', line)):
					pt['project'] = line.strip().split(':')[1]
				if(re.search('Local time:', line)):
					pt['dateTime'] = line.strip().split(':', 1)[1]
				if(re.search('epochTime:', line)):
					pt['epochTime'] = float(line.strip().split(':', 1)[1])
				if(re.search('startingwfc:', line)):
					pt['startingwfc'] = line.strip().split(':', 1)[1]
				if(re.search('startingpot:', line)):
					pt['startingpot'] = line.strip().split(':', 1)[1]
				if(re.search('the Fermi energy is', line)):
					m = re.search('the Fermi energy is[ ]*(?P<fermiLevel>[0-9.]+)[ ]*ev', line)
					d = m.groupdict()
					pt['fermiLevel'] = float(d['fermiLevel'])
				if(re.search(r'highest occupied, lowest unoccupied level \(ev\):', line)):
					m = re.search(r'highest occupied, lowest unoccupied level \(ev\):[ ]*(?P<vbm>[0-9.]+)[ ]+(?P<cbm>[0-9.]+)', line)
					d = m.groupdict()
					pt['vbm'] = float(d['vbm'])
					pt['cbm'] = float(d['cbm'])
					pt['bandGap'] = pt['cbm'] - pt['vbm']
				if(re.search(r'highest occupied level \(ev\):', line)):
					m = re.search(r'highest occupied level \(ev\):[ ]*(?P<vbm>[0-9.]+)', line)
					d = m.groupdict()
					pt['vbm'] = float(d['vbm'])
				if(re.search(r'unit-cell volume', line)):
					m = re.search(r'unit-cell volume[ ]*=[ ]*(?P<volume>[0-9.]+)[ ]+', line)
					d = m.groupdict()
					pt['V'] = float(d['volume'])*0.001*(self.bohr2ang**3)
				if(re.search('Exchange-correlation[ ]*=', line)):
					m = re.search('Exchange-correlation[ ]*=[ ]*(?P<XCfunc>[a-zA-Z0-9-]+)[ ]+', line)
					d = m.groupdict()
					pt['XCfunc'] = d['XCfunc']
				if(re.search('number of electrons[ ]*=', line)):
					m = re.search('number of electrons[ ]*=[ ]*(?P<electrons>[0-9.]+)', line)
					d = m.groupdict()
					pt['valenceElectrons'] = float(d['electrons'])
					self.valenceElectrons = pt['valenceElectrons']
			if log:
				pt['calcType'] = log
			pt['nat'] = self.nat
			self.data.append(pt)
		if log:
			if os.path.isfile(self.datalog):
				with open(self.datalog, 'r') as file:
					dataList = json.loads(file.read())
			else:
				dataList = []
			dataList = dataList + self.data
			with open(self.datalog, 'w+') as file:
				file.write(json.dumps(dataList))
		return self.data
	def write(self, calcType):
		self.clearDir(self.inputFileFolder, False)
		if self.ibrav == 0 and len(self.crystal.translations) == 0:
			raise Exception("You must manually enter cell tranlation vectors (cellParams) when using the free lattice (ibrav = 0). Set cellParams to a 3 element list of three element cartesian translation vectors (floats). ")
		inputFile = open(os.path.join(self.ipFileDir, '%s-%s.in' % (self.postfix, calcType)),'w+')
		inputFile.truncate()
		inputFile.write(" &CONTROL\n")
		inputFile.write("  calculation = '%s',\n" % self.calculation)
		if self.restart:
			inputFile.write("  restart_mode = 'restart',\n")
		inputFile.write("  outdir = '%s',\n" % self.outdir)
		inputFile.write("  pseudo_dir = '%s',\n" % self.pseudoDir)
		inputFile.write("  prefix = '%s',\n" % self.prefix)
		inputFile.write("  verbosity = '%s',\n" % self.verbosity)
		if self.wfCollect:
			inputFile.write("  wf_collect = .TRUE.,\n")
		if self.nstep:
			inputFile.write("  nstep = %i,\n" % self.nstep)
		inputFile.write(" /\n")
		inputFile.write(" &SYSTEM\n")
		inputFile.write("  ibrav = %i,\n" % self.ibrav)
		inputFile.write("  celldm(1) = %f,\n" % (self.crystal.a/self.bohr2ang))
		inputFile.write("  nat = %i,\n" % self.nat)
		inputFile.write("  ntyp = %i,\n" % self.ntyp)
		inputFile.write("  ecutwfc = %d,\n" % self.ecutwfc)
		inputFile.write("  ecutrho = %d,\n" % self.ecutrho)
		if self.hubbardU:
			inputFile.write("  lda_plus_u = .TRUE.,\n")
			if self.ldaUtype1 or self.lspinorb:
				inputFile.write("  lda_plus_u_kind = 1,\n")
			for key, value in iteritems(self.hubbardU):
				for index, element in enumerate(self.crystal.elements):
					if element in self.hubbardU:
						inputFile.write("  Hubbard_U(%i) = %f,\n" % (index+1, value))
		if self.noSym:
			inputFile.write("  nosym = .TRUE.,")
		if self.nbnd:
			inputFile.write("  nbnd = %i, \n" % (self.nbnd*self.bandMultiplier))
		if self.occupations:
			inputFile.write("  occupations = '%s' ,\n" % self.occupations)
		if self.degauss:
			inputFile.write("  degauss = %f, \n" % self.degauss)
		if self.input_dft:
			inputFile.write("  input_dft = '%s', \n" % self.input_dft)
		if self.charge:
			inputFile.write("  tot_charge = %f, \n" % self.charge)
		if self.lspinorb:
			inputFile.write("  lspinorb = .TRUE., \n")
			inputFile.write("  noncolin = .TRUE., \n")
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
				if self.cellDoFree:
					inputFile.write(" cell_dofree = '%s'\n" % self.cellDoFree)
				inputFile.write(" /\n")
		inputFile.write("ATOMIC_SPECIES\n")
		for symbol,element in iteritems(self.crystal.elements):
			if symbol not in self.pseudoFiles:
				exit('Missing pseudopotential file for %s' % symbol)
			inputFile.write(" %s %f %s\n" % (symbol, element.mass, self.pseudoFiles[symbol]))
		inputFile.write("ATOMIC_POSITIONS alat\n")
		for atom in self.crystal.basis:
			inputFile.write(" %s %f %f %f %i %i %i\n" % (atom[0], atom[1][0], atom[1][1], atom[1][2], atom[2][0], atom[2][1], atom[2][2]))
		if self.kPath:
			inputFile.write("K_POINTS crystal_b\n")
			inputFile.write(" %i\n" % len(self.kPath))
			for pt in self.kPath:
				inputFile.write("  %f %f %f %i \n" % (pt[0], pt[1], pt[2], self.kPathDensity))
		else:
			inputFile.write("K_POINTS automatic\n")
			shift = '1 1 1' if self.shiftK else '0 0 0'
			inputFile.write("%d %d %d %s\n" % (self.kxyz[0], self.kxyz[1], self.kxyz[2], shift))
		if self.ibrav == 0:
			inputFile.write("CELL_PARAMETERS alat\n")
			for v in self.crystal.translations:
				inputFile.write(" %f %f %f\n" % (v[0], v[1], v[2]))
		if (self.calculation == 'relax' or self.calculation == 'vc-relax') and self.constraints:
			inputFile.write("CONSTRAINTS\n")
			if self.constraints == 'zb-bennett':
				inputFile.write("1\n")
				inputFile.write("'bennett_proj' 0 1.0 1.0 1.0 \n")
		inputFile.close()
	def writeDOS(self):
		inputFile = open(self.ipDOS,'w+')
		inputFile.truncate()
		#inputFile.write(" &DOS\n")
		inputFile.write(" &PROJWFC\n")
		inputFile.write("  outdir = '%s',\n" % self.outdir)
		inputFile.write("  prefix = '%s',\n" % self.prefix)
		inputFile.write("  DeltaE = %f,\n" % self.dosDeltaE) #DeltaE
		if self.dosEmin:
			inputFile.write("  Emin = %f,\n" % self.dosEmin)
		if self.dosEmax:
			inputFile.write("  Emax = %f,\n" % self.dosEmax)
		inputFile.write("  filpdos = '%s',\n" % self.fildos)
		inputFile.write(" /\n")
	def writeBands(self):
		inputFile = open(self.ipBands,'w+')
		inputFile.truncate()
		inputFile.write(" &BANDS\n")
		inputFile.write("  outdir = '%s',\n" % self.outdir)
		inputFile.write("  prefix = '%s',\n" % self.prefix)
		inputFile.write("  filband = '%s',\n" % self.filband)
		if self.lspinorb:
			inputFile.write("  spin_component = 1,\n")
		inputFile.write(" /\n")
	def run(self, cmdStr, filepath):
		"""Run a Quantum Espresso routine, reading stdout for errors before writing to file"""
		result = GenericClass()
		result.jobDone = False
		result.errorEncountered = False
		result.errorType = None 
		with open(filepath, 'wb+') as opFile:
			#shellCommand = "%s %s > %s" % (os.path.join(self.pgmDir, cmdStr), self.ipBands, os.path.join(self.bandPath,'%s-bandsx.out' % self.postfix))
			with Popen(cmdStr, stdout=PIPE) as p:
				printing = False
				for line in p.stdout:
					if re.search(b'%%%%%%%%%%%%%%%', line):
						if result.errorEncountered:
							printing = False
						else:
							result.errorEncountered = True
							self.log(self.highlight('ERROR ENCOUNTERED'))
							printing = True
					if re.search(b'Error in routine d_matrix', line):
						result.errorType = self.lastError = 'dMatrix'
					if re.search(b'use lda_plus_u_kind = 1', line):
						if self.salvageNumber < 2:
							self.log(self.highlight('Wrong LDA+U type error. Attempting to restart with corrected type.'))
							self.ldaUtype1 = True
							return self.run(cmdStr, filepath)
						else:
							result.errorEncountered = True
							result.errorType = 'ldaUtype'
					if re.search(b'JOB DONE.', line):
						result.jobDone = True
					if printing:
						self.log(line.decode('utf-8'))
					opFile.write(line)
		if not result.jobDone:
			if os.path.isfile('CRASH'):
				result.errorEncountered = True
				with open('CRASH', 'r') as f:
					for line in f.readlines():
						if re.search('lone vector', line):
							result.errorType = self.lastError = 'loneVector'
		return result
	def runPw(self, opFilename, calcType):
		"""Runs pw.x with the current input file and writes the output to opFilename. self.mpi is a 5-element array that reads [np,ni,nk,nb,nt,nd] (see QE user guide page 23). if self.mpi is set, the calculation will be run in parallel using mpirun"""
		pickleFile = open(os.path.join(self.workingDir, '%s-io.pkl' % self.postfix),'wb')
		pickleFile.truncate()
		pickle.dump(self, pickleFile, pickle.HIGHEST_PROTOCOL)
		if config.system.os == 'Windows':
			pgm = 'pw'
		if config.system.os == 'Linux':
			pgm = 'pw.x'		
		if self.mpi:
			nprocs = multiprocessing.cpu_count()
			if len(self.mpi) < 6:
				print('Your mpi parameters array is too short')
				exit()
			#if self.mpi[0] < nprocs:
				#print('Warning: you are running in parallel, but not utilizing all of your processors')
			if self.mpi[0] > nprocs:
				print('You have specified more processors than are available')
				exit()
			if self.mpi[1]*self.mpi[2] > self.mpi[0]:
				print('Your mpi paramters specify too many images/pools. Use constraint ni*nk <= np')
				exit()
			if self.mpi[5] > self.mpi[0]/(self.mpi[1]*self.mpi[2]):
				print('You have specified too many processors for your diagonalization group.')
				exit()
			if not self.mpi[5] == int(math.sqrt(self.mpi[5]))**2:
				print('The diagonalization group size must be a perfect square')
				exit()
			#if not isinstance(self.mpi[2]/self.mpi[4], int):
				#print('Warning: your mpi params are not optimized. Try constraint nk/nt = Int')
			#cmdStr = ' '.join(['mpirun', '-np', str(self.mpi[0]), os.path.join(self.pgmDir,'pw.x'), '-ni', str(self.mpi[1]), '-nk', str(self.mpi[2]), '-nb', str(self.mpi[3]), '-nt', str(self.mpi[4]), '-nd', str(self.mpi[5]), '-i', os.path.join(self.ipFileDir, '%s-%s.in' % (self.postfix, calcType)), '>', opFilename])
			cmdStr = ['mpirun', '-np', str(self.mpi[0]), os.path.join(self.pgmDir,pgm), '-ni', str(self.mpi[1]), '-nk', str(self.mpi[2]), '-nb', str(self.mpi[3]), '-nt', str(self.mpi[4]), '-nd', str(self.mpi[5]), '-i', os.path.join(self.ipFileDir, '%s-%s.in' % (self.postfix, calcType))]
		else:
			#cmdStr = "%s %s > %s" % (self.pgmDir+'/'+pgm, os.path.join(self.ipFileDir, '%s-%s.in' % (self.postfix, calcType)), opFilename)
			cmdStr = [os.path.join(self.pgmDir, pgm), '-i', os.path.join(self.ipFileDir, '%s-%s.in' % (self.postfix, calcType))]
		result = self.run(cmdStr, opFilename)
		self.appendData(opFilename)
		return result
	def runDOS(self, Ecenter=False, Espread=False, Emin=False, Emax=False, deltaE = 0.01, plot=False):
		if Ecenter and Espread:
			self.dosEmin = Ecenter - Espread
			self.dosEmax = Ecenter + Espread
			self.dosEcenter = Ecenter
		elif Ecenter or Espread:
			print('Ecenter and Espread must be used together')
			exit()
		else:
			self.dosEmin = Emin
			self.dosEmax = Emax
		self.dosDeltaE = deltaE
		self.writeDOS()
		self.clearDir(self.dosFileDir)
		if config.system.os == 'Windows':
			pgm = 'projwfc'
		if config.system.os == 'Linux':
			pgm = 'projwfc.x'
		self.log('projwfc.x calculation commenced at %s (server time)' % time.asctime())
		cmdStr = cmdStr = [os.path.join(self.pgmDir,pgm), '-i', self.ipDOS]
		if self.mpi:
			cmdStr = ['mpirun', '-np', str(self.mpi[0])] + cmdStr
		result = self.run(cmdStr, os.path.join(self.dosPath, '%s-dos.out' % self.postfix))
		if result.errorEncountered:
			return result
		self.plotDOS(plot = plot)
		return result
	def runBands(self, plot=False, bandLabels = False):
		""" Create band output file for plotting on the current object using last pw output as data for bands.x"""
		self.writeBands()
		self.clearDir(self.bandFileDir)
		if config.system.os == 'Windows':
			pgm = 'bands'
		if config.system.os == 'Linux':
			pgm = 'bands.x'
		cmdStr = [os.path.join(self.pgmDir,pgm), '-i', self.ipBands]
		if self.mpi:
			cmdStr = ['mpirun', '-np', str(self.mpi[0])] + cmdStr
		result = self.run(cmdStr, os.path.join(self.bandPath,'%s-bandsx.out' % self.postfix))
		self.plotBands(bandLabels = bandLabels, plot=plot)
		return result
	def plotBands(self, bandLabels = False, plot=False):
		"""Plot the band file created by runBands"""
		data = self.parseResults(self.bandPwFileDir)
		if isinstance(data[-1]['vbm'], float):
			self.vbm = data[-1]['vbm']
		offset = self.vbm if self.vbm else 0
		with open(self.filband, 'r') as f:
			lines = f.readlines()
			m = re.search('&plot[ ]*nbnd=[ ]*(?P<nbnd>[0-9]+),[ ]*nks=[ ]*(?P<nks>[0-9]+)', lines[0])
			d = m.groupdict()
			nbnd = int(d['nbnd'])
			dataLines = int(math.ceil(nbnd/10))
			nks = int(d['nks'])
			Y = [[] for x in range(nbnd)]
			X = []
			dataFileArray = [' '.join(['i x y z',' '.join([str(x) for x in range(nbnd)])])]
			for i, line in enumerate(lines[1:-1:dataLines + 1]):
				data = []
				for j in range(1, dataLines+1):
					data.append(lines[i*(dataLines + 1) + j + 1].strip())
				data = ' '.join(data)
				a = data.split()
				for index, value in enumerate(a):
					Y[index].append(float(value) - offset)
				dataFileArray.append(' '.join(['%i %s' % (i, line.strip()), data]))
				X.append(i)
		with open(os.path.join(self.bandPath, '%s-bands.data' % self.postfix), 'w+') as f:
			f.write('\n'.join(dataFileArray))
		fig = plt.figure()
		fig.patch.set_facecolor('white')	
		#plt.legend(loc="upper center", fancybox=True, ncol=1, fontsize=13, numpoints = 1)
		plt.suptitle('L-$\Gamma$-X Bands - %s' % self.postfix, fontsize=20)
		#plt.xlabel(r'Lattice constant ($\AA$)', fontsize=18)
		plt.ylabel('Energy (eV)', fontsize=18)
		if bandLabels:
			labelx = []
			labels = []
			for index, label in enumerate(bandLabels):
				labelx.append(index*self.kPathDensity)
				labels.append(label)
			ax = plt.axes()
			ax.set_xlim(0,(len(self.kPath)-1)*self.kPathDensity)
			plt.xticks(labelx,labels)
		for band in Y:
			plt.plot(X,band)
		fig.savefig(os.path.join(self.workingDir, '%s-bands.png' % self.postfix),format='png')
		if config.system.hasDisplay and plot:
			plt.ion()
			plt.show()
			if self.waitAtImages:
				raw_input('Press Enter to Continue')
		else:
			plt.close()
	def plotDOS(self, imgName=False, plot=False):
		""" Create a plot from the data files of runDOS output"""
		data = self.parseResults(self.dosPwFileDir)
		if isinstance(data[-1]['vbm'], float):
			self.vbm = data[-1]['vbm']
		if isinstance(data[-1]['fermiLevel'], float):
			self.fermiLevel = data[-1]['fermiLevel']
		offset = self.vbm if self.vbm else 0
		if not self.vbm:
			offset = self.fermiLevel if self.fermiLevel else 0
		data = {}
		for filename in os.listdir(self.dosPath):
			if re.search('.pdos_tot$', filename):
				key = 'total'
			elif re.search('.pdos_atm', filename):
				m = re.search(r'.pdos_atm#[0-9]+\((?P<symbol>[a-zA-Z]+)\)_wfc#[0-9]+\((?P<orbital>[a-z_.0-9]+)\)', filename)
				d = m.groupdict()
				key = '%s-%s' % (d['symbol'],d['orbital'])
			else:
				continue
			data[key] = {'x':[],'y':[]}
			with open(os.path.join(self.dosPath, filename), 'r') as f:
				lines = f.readlines()
				for line in lines[1:-1]:
					a = line.strip().split()
					data[key]['x'].append(float(a[0]) - offset)
					data[key]['y'].append(float(a[1]))
		fig = plt.figure()
		fig.patch.set_facecolor('white')
		for key, d in iteritems(data):
			s = plt.plot(d['x'], d['y'])
			plt.setp(s, label = key)
		plt.legend(loc="upper right", fancybox=True, ncol=1, fontsize=13, numpoints = 1)
		plt.suptitle('Density of States', fontsize=20)
		plt.xlabel(r'Energy (eV)', fontsize=18)
		plt.ylabel('DOS', fontsize=18)
		fig.savefig(os.path.join(self.workingDir, '%s-dos.png' % self.postfix),format='png')
		if config.system.hasDisplay and plot:
			plt.ion()
			plt.show()
			if self.waitAtImages:
				raw_input('Press Enter to Continue')
		else:
			plt.close()
		outData = [[]]
		outData[0] = ['energy'] + [str(x) for x in data['total']['x']]
		for key, d in iteritems(data):
			outData.append([key] + [str(y) for y in d['y']])
		with open(os.path.join(self.dosPath, '%s-dos.data' % self.postfix), 'w+') as f:
			f.write('\n'.join([' '.join(z) for z in zip(*outData)]))
	def fitLattice(self, plot=False):
		""" Plots the output of self.parseResults after a lattice optimization """
		self.parseResults(self.latticeFileDir)
		x = []
		y = []
		pts = sorted(self.data, key=lambda k: float(k['latticeConstant']))
		if len(pts) < 5:
			print("You don't appear to have enough data points. Be sure to collect at least 5 data points and to call qeio.parseResults before fitting the lattice data.")
			return False
		for pt in pts:
			y.append(pt['energy'])
			x.append(pt['latticeConstant'])
		dataFileArray = ['lattice-constant(Angstrom) energy(eV)'] + ['%s %s' % (j,k) for j,k in zip(x,y)]
		with open(os.path.join(self.sub2path(self.latticeFileDir), '%s-lattice.data' % self.postfix), 'w+') as f:
			f.write('\n'.join(dataFileArray))
		def quad(a,E0,B,a0):
			return E0 + B*(a - a0)**2
		try:
			quadParams, pcov = curve_fit(quad, x, y, p0=[y[1],1,x[1]])
		except RuntimeError:
			self.lastError = 'quadFit'
			return False
		def BMeos(a,a0,dB,E0,B0):
			A = (a0/a)**2.
			return E0 + (9./16.)*a0**3.*la.det([self.crystal.translations[0], self.crystal.translations[1], self.crystal.translations[2]])*B0*((A - 1.)**3.*dB + (A - 1.)**2*(6. - 4.*A))
		try:
			BMparams, pcov = curve_fit(BMeos, x, y, p0=[x[1],4.,y[1],0.2])
		except RuntimeError:
			self.lastError = 'BMfit'
			return False
		B0 = 160.518*BMparams[3]
		a0 = BMparams[1]
		denseX = []
		tmpX = min(x)
		deltaX = (max(x) - min(x))/30
		while tmpX < max(x):
			denseX.append(tmpX)
			tmpX += deltaX
		fig = plt.figure()
		fig.patch.set_facecolor('white')			
		dataMarkers = plt.plot(x,y)
		aveX = (min(x) + max(x))/2
		aveY = (min(y) + max(y))/2
		ySpread = (max(y) - min(y))
		plt.setp(dataMarkers, color='k', marker='s', linestyle = 'none', markersize = 8.0, label = 'Data')
		quadY = [quad(a, quadParams[0], quadParams[1], quadParams[2]) for a in denseX]
		quadLine = plt.plot(denseX, quadY)
		plt.setp(quadLine,color='k', marker=None, linestyle='solid', linewidth=1.0, label = 'Quadratic')
		BMy = [BMeos(a, BMparams[0], BMparams[1], BMparams[2], BMparams[3]) for a in denseX]
		BMline = plt.plot(denseX, BMy)
		plt.setp(BMline,color='k', marker=None, linestyle='dotted', linewidth=1.0, label = 'BM')
		plt.legend(loc="upper center", fancybox=True, ncol=1, fontsize=13, numpoints = 1)
		plt.suptitle('DFT energy v. Lattice Constant for FCC CdTe', fontsize=20)
		plt.xlabel(r'Lattice constant ($\AA$)', fontsize=18)
		plt.ylabel('DFT energy (Ry)', fontsize=18)
		plt.subplots_adjust(left = 0.2, right = 0.9, top = 0.85, bottom = 0.2)
		plt.text(aveX, aveY - ySpread*0.05, 'B = %.2f GPa' % B0)
		plt.text(aveX, aveY + ySpread*0.05, 'a = %.3f $\AA$' % a0)
		filePath = '.'.join([os.path.join(self.workingDir, self.postfix), '-lattice.png'])
		fig.savefig(filePath,format='png')
		if config.system.hasDisplay and plot:
			plt.ion()
			plt.show()
		else:
			plt.close()
		returnDict = {}
		returnDict['BMfit'] = a0
		returnDict['bulkMod'] = B0
		returnDict['quadFit'] = quadParams[2]
		return returnDict
	def fitCutoff(self, qeParam, subDir, plot=False):
		""" Plots and fits trendlines to data from energy cutoff convergence calculation"""
		self.parseResults(subDir)
		x = []
		y = []
		t = []
		pts = sorted(self.data, key=lambda k: float(k[qeParam]))
		for pt in pts:
			y.append(pt['energy'])
			x.append(pt[qeParam])
			t.append(pt['wallTime'])
		return self.fitConvergence(x, y, t, subDir, xLabel=qeParam+' Cutoff', plot=plot)
	def fitKgrid(self, subDir, plot=False):
		""" Plots and fits trendlines to data from density energy cutoff convergence calculation"""
		self.parseResults(subDir)
		x = []
		y = []
		t = []
		xTickLabels = []
		kGridX = []
		pts = sorted(self.data, key=lambda k: float(k['kPts']))
		for pt in pts:
			y.append(pt['energy'])
			x.append(float(pt['kGridDict']['x']) + float(pt['kGridDict']['y']) + float(pt['kGridDict']['z']))
			xTickLabels.append('[%s,%s,%s]' % (pt['kGridDict']['x'], pt['kGridDict']['y'], pt['kGridDict']['z']))
			kGridX.append([int(pt['kGridDict']['x']), int(pt['kGridDict']['y']), int(pt['kGridDict']['z'])])
			t.append(pt['wallTime'])
		return self.fitConvergence(x, y, t, subDir, xLabel = 'k-Grid Divisions', xTickLabels=xTickLabels, plot=plot, kGridX = kGridX)
	def fitConvergence(self, x, y, t, subDir, xLabel='', xTickLabels=False, plot=False, kGridX = False):
		fileString = 'kgrid' if subDir == self.kGridFileDir else 'ecutwfc' if subDir == self.ecutwfcFileDir else 'ecutrho'
		fileX = x
		headerX = fileString
		if fileString == 'kgrid':
			fileX = ['%i %i %i' % (i,j,k) for i,j,k in kGridX]
			headerX = 'Mx My Mz'
		with open(os.path.join(self.sub2path(subDir), '%s-%s.data' % (self.postfix, fileString)), 'w+') as f:
			f.write('\n'.join(['%s energy(eV) time(s)' % headerX] + ['%s %f %f' % (str(i),j,k) for i,j,k in zip(fileX,y,t)]))
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
		filePath = '.'.join([os.path.join(self.workingDir, '%s-%s' % (self.postfix,fileString)), 'png'])
		fig.savefig(filePath,format='png')
		returnDict = {}
		returnDict['files'] = []
		returnDict['files'].append(filePath)
		if config.system.hasDisplay and plot:
			plt.ion()
			plt.show()
			if self.waitAtImages:
				raw_input('Press Enter button to continue')
		else:
			plt.close()
		return returnDict
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
	def mkDir(self, subDir):
		""" Create empty directory, if it doesn't exist. Same as clearDir() with sweep=False"""
		dirPath = self.sub2path(subDir)
		if not os.path.isdir(dirPath):
			try:
				os.mkdir(dirPath)
			except:
				print('Error creating output directory\n')
				print('Error data:'+str(sys.exc_info()))
				return False
	def rmDir(self, subDir):
		"""Delete directory and all of it's contents"""
		dirPath = self.sub2path(subDir)
		if os.path.isdir(dirPath):
			shutil.rmtree(dirPath)
	def clearDir(self, subDir, sweep=True):
		""" Clears directory. Creates empty directory if it doesn't exist"""
		dirPath = self.sub2path(subDir)
		if not os.path.isdir(dirPath):
			try:
				os.mkdir(dirPath)
			except:
				print('Error creating output directory\n')
				print('Error data:'+str(sys.exc_info()))
				return False
		if sweep:
			for file in os.listdir(dirPath):
				filePath = os.path.join(dirPath, file)
				try:
					if os.path.isfile(filePath):
						os.unlink(filePath)
				except:
					print('Error clearing directory')
					print('Error data:'+sys.exc_info())
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
			print('Kill switch activated. If any calculations are running for the provided project, they should end shortly.')
		else:
			print('There does not seem to be any project with that name.')
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
			print("There doesn't appear to be a calculation to restart")
	def log(self, entry):
		if not entry:
			return
		print(entry)
		with open(self.logPath, 'a') as logFile:
			logFile.write(str(entry)+'\n')
	def logValue(self, values):
		""" Log a value in a .dictionary file, used for holding important values. """
		data = OrderedDict()
		if os.path.isfile(self.dictionaryPath):
			with open(self.dictionaryPath, 'r') as f:
				for line in f.readlines():
					k, v = line.strip().split()
					data[k] = v
		data['project'] = self.project
		data['postfix'] = self.postfix
		for k, v in iteritems(values):
			data[k] = v
		with open(self.dictionaryPath, 'w+') as f:
			for k, v in iteritems(data):
				f.write('%s %s\n' % (k,v))
	def readDictionary(self):
		"""Look for a .dictionary,and if it exists, create a python dictionary type object"""
		valueDict = {}
		if os.path.isfile(self.dictionaryPath):
			with open(self.dictionaryPath, 'r') as f:
				for line in f.readlines():
					k, v = line.split()
					valueDict[re.sub(r'\([a-zA-Z0-9]+\)', '', k)] = v
		return valueDict
	def highlight(self, text):
		return "\033[93m%s\033[0m" % text
	def report(self, message, email=False, subject = '', text=False, files=False, image=False):
		if email:
			send.email(subject, message, files=files)
		if text:
			send.text(message)
	def exportData(self):
		if os.path.isfile(self.datalog):
			with open(self.datalog, 'r') as file:
				dataList = json.loads(file.read())
			#print(str(dataList))
			with open(self.csvPath, 'w+') as file:
				writer = csv.DictWriter(file, self.flattenData(dataList[0]).keys())
				writer.writeheader()
				for dataPt in dataList:
					writer.writerow(self.flattenData(dataPt))
		else:
			print('No data to export.')
	def flattenData(self, dataPt):
		opDict = {}
		if 'mpi' not in dataPt:
			dataPt['mpi'] = {'np':'1','ni':'1','nk':'1','nb':'1','nt':'1','nd':'1'}
		for key in dataPt:
			if isinstance(dataPt[key], dict):
				for subkey in dataPt[key]:
					opDict['.'.join([key,subkey])] = dataPt[key][subkey]
			else:
				opDict[key] = dataPt[key]			
		return opDict
	def initDirs(self):
		for path in [self.workingDir]:
			if not os.path.isdir(path):
				try:
					os.mkdir(path)
				except:
					print('Error creating directory: %s' % path)
					print('Error data:'+sys.exc_info())
			
			
