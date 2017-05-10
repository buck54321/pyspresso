from compatibility import *
import functools
from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *
from PyQt4 import QtGui, QtCore
from PyQt4.QtOpenGL import *
import geometry as geo
import numpy.linalg as la
import numpy as np
import random


class MainWindow(QtGui.QWidget):
	def __init__(self):
		super(MainWindow, self).__init__()
		self.controlBox = ControlWindow(self)
		self.glBox = glWidget(self)
		self.mainLayout = QtGui.QHBoxLayout()
		self.mainLayout.addWidget(self.glBox)
		self.mainLayout.addWidget(self.controlBox)
		self.setLayout(self.mainLayout)
		self.controlBox.bind()
	def maximize(self):
		self.showMaximized()
		self.relay()
	def relay(self):
		self.glBox.setGeometry(0,0,self.width() - 200, self.height())
		self.controlBox.setGeometry(self.width() - 200,0,200,self.height())
	def resizeEvent(self,e):
		self.glBox.setGeometry(0,0,self.width() - 200, self.height())
		self.controlBox.setGeometry(self.width() - 200,0,200,self.height())
	def loadAtoms(self, atoms, translations):
		self.glBox.loadAtoms(atoms, translations)

class ControlWindow(QtGui.QWidget):
	def __init__(self, parent):
		super(QtGui.QWidget, self).__init__()
		self.parent = parent
		self.glBox = False
		self.setMinimumSize(200,200)
		self.setFixedWidth(200)
		self.layout = QtGui.QGridLayout()
		self.layout.setAlignment(QtCore.Qt.AlignTop)
		self.xPosLabel, self.xPosValue = self.addValue('x position','0')
		self.layout.addWidget(self.xPosLabel, 0, 0)
		self.layout.addWidget(self.xPosValue, 0, 1)
		self.yPosLabel, self.yPosValue = self.addValue('y position', '0')
		self.layout.addWidget(self.yPosLabel, 1, 0)
		self.layout.addWidget(self.yPosValue, 1, 1)
		self.randomizeColorButton = QtGui.QPushButton('Randomize Colors', self)
		self.randomizeColorButton.setFixedWidth(150)
		self.randomizeColorButton.setFixedHeight(30)
		self.layout.addWidget(self.randomizeColorButton, 2, 0, 1, 2)
		self.edgeAtomButton = QtGui.QPushButton('Show Edge Atoms', self)
		self.edgeAtomButton.setFixedWidth(150)
		self.edgeAtomButton.setFixedHeight(30)
		self.edgeAtomButton.setCheckable(True)
		self.layout.addWidget(self.edgeAtomButton, 3, 0, 1, 2)

		self.setLayout(self.layout)

	def bind(self):
		self.glBox = self.parent.glBox
		self.randomizeColorButton.clicked.connect(self.glBox.randomizeColors)
		self.edgeAtomButton.clicked.connect(functools.partial(self.glBox.toggleEdgeAtoms, self.edgeAtomButton.isChecked())) # isChecked is always returning false for some reason

	def addValue(self, lab, val='0'):
		label = QtGui.QLabel(lab, self)
		label.setStyleSheet('font-size:20px;')
		label.setFixedWidth(100)
		label.setFixedHeight(30)
		label.setAlignment(QtCore.Qt.AlignRight)
		value = QtGui.QLabel(val, self)
		value.setStyleSheet('font-size:15px;background-color:white;')
		value.setFixedWidth(50)
		value.setFixedHeight(30)
		value.setAlignment(QtCore.Qt.AlignCenter)
		return (label,value)


class glWidget(QGLWidget):
	def __init__(self, parent):
		QGLWidget.__init__(self, parent)
		self.parent = parent
		self.w = self.parent.size().width() - 200
		self.h = self.parent.size().height()
		self.aspectRatio = self.w/self.h
		self.rotFactorX = self.h/360*np.pi/180
		self.rotFactorY = self.w/360*np.pi/180
		self.includeEdgeRepeats = True
		self.setMinimumSize(200, 200)
		self.setGeometry(0, 0, self.w, self.h)
		self.atoms = []
		self.edgeRepeats = []
		self.translations = []
		self.centerV = [0.,0.,0.]
		self.cellSizeFactor = 0.7
		self.atomSizeFactor = 0.07
		self.atomSize = self.atomSizeFactor*self.cellSizeFactor
		self.quadric = gluNewQuadric()
		self.mouseIsDown = False
		self.setMouseTracking(True)
		self.xPosDisplay = self.parent.controlBox.xPosValue
		self.yPosDisplay = self.parent.controlBox.yPosValue
		self.xDown = 0
		self.yDown = 0
		self.trackingX = np.array([1.,0.,0.])
		self.trackingY = np.array([0.,1.,0.])
		self.yRotation = 0
		self.xRotation = 0
		self.Mrot = geo.I
		self.MinvList = self.a32l16(geo.I)
		self.atomColors = {}
		self.xyz = []
		self.cellMaxDim = 1
		self.cylinderQuadric =gluNewQuadric() # // create a new Quadric object
		gluQuadricDrawStyle(self.cylinderQuadric, GLU_FILL)  # // set the drawing style
		gluQuadricNormals(self.cylinderQuadric, GLU_SMOOTH)  # // create normals for shading   
		self.cylinderRadius = self.atomSize/4


	def paintGL(self):
		glMatrixMode(GL_MODELVIEW)
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
		glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, (1.,1.,1.,1.))
		glPolygonMode(GL_FRONT, GL_FILL);

		for atom in self.atoms + self.edgeRepeats:
			glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, self.atomColors[atom[0]])
			glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, self.atomColors[atom[0]])
			glLoadIdentity()
			glMultMatrixf(self.MinvList)
			glTranslatef(atom[1][0]*self.cellSizeFactor, atom[1][1]*self.cellSizeFactor, atom[1][2]*self.cellSizeFactor)
			gluSphere(self.quadric,self.atomSize,32,32)
		
		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, (0.,0.25,0.25,1.))
		glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, (0.,0.25,0.25,1.))
		for i, t in enumerate(self.translations):
			others = [x for x in self.translations]
			others.pop(i)
			others = [v for l, m, v in others]
			others.append(np.array(others[0])+np.array(others[1]))
			others.append([0.,0.,0.])
			for other in others:
				glLoadIdentity()
				glMultMatrixf(self.MinvList)
				glTranslatef(*np.array(other)*self.cellSizeFactor)
				glTranslatef(*self.centerV*self.cellSizeFactor*-1)
				glMultMatrixf(t[1])
				gluCylinder(self.cylinderQuadric, self.cylinderRadius, self.cylinderRadius, t[0]*self.cellSizeFactor, 16, 5)


		glFlush()

	def initializeGL(self):
		glClearColor(1.0,1.0,1.0,0.)
		glClearDepth(1.0)              
		glDepthFunc(GL_LESS)
		glLineWidth(5.0)
		glEnable(GL_DEPTH_TEST)
		glShadeModel(GL_SMOOTH)
		
		glEnable(GL_LIGHTING)
		glEnable(GL_LIGHT0)
		glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE)
		glLightfv(GL_LIGHT0, GL_POSITION, (50.,50.,100.,1.))
		glLightfv(GL_LIGHT0, GL_DIFFUSE, (1.,1.,1.,1.))
		glLightfv(GL_LIGHT0, GL_SPECULAR, (0.2,0.2,0.2,1.))
		
		self.setScene()

	def resizeGL(self, w, h):
		self.setScene()

	def setScene(self):
		self.w = self.parent.size().width() - 200
		self.h = self.parent.size().height()
		self.aspectRatio = self.w/self.h
		self.rotFactorX = 360/self.h*np.pi/180
		self.rotFactorY = 360/self.w*np.pi/180
		glMatrixMode(GL_PROJECTION)
		glLoadIdentity()
		glViewport(0,0,self.w,self.h)
		glLoadIdentity()
		glOrtho(-0.5, 0.5, -0.5/self.aspectRatio, 0.5/self.aspectRatio, -500, 500)
		glMatrixMode(GL_MODELVIEW)
		glLoadIdentity()
		gluLookAt(0.0, 0.0, 1.0, 0.0, 0.0, -1.0, 0.0, 1.0, 0.0) 

	def mousePressEvent(self, e):
		self.mouseIsDown = True
		self.xDown = e.x()
		self.yDown = e.y()
	def mouseReleaseEvent(self, e):
		self.mouseIsDown = False
	def mouseMoveEvent(self, e):
		if self.mouseIsDown:
			y = e.y()
			x = e.x()
			ty = (x - self.xDown)*self.rotFactorY
			tx = (y - self.yDown)*self.rotFactorX
			self.trackingX = geo.rotOnY(self.trackingX, ty)
			self.trackingX = geo.rotOnX(self.trackingX, tx)
			self.trackingY = geo.rotOnY(self.trackingY, ty)
			self.trackingY = geo.rotOnX(self.trackingY, tx)
			vx = self.trackingX
			vy = self.trackingY
			self.Mrot = rx = geo.Rx(-1*geo.yzy(vx))
			vx = np.inner(rx, vx)
			vy = np.inner(rx, vy)
			rz = geo.Rz(-1*geo.xyx(vx))
			self.Mrot = np.matmul(rz,self.Mrot)
			vx = np.inner(rz, vx)
			vy = np.inner(rz, vy)
			rx = geo.Rx(-1*geo.yzy(vy))
			self.Mrot = np.matmul(rx, self.Mrot)
			self.MinvList = self.a32l16(la.inv(self.Mrot))
			self.xDown = x
			self.yDown = y
			self.updateGL()
			#self.xPosDisplay.setText('%i' % x)
			#self.yPosDisplay.setText('%i' % y)
	def loadAtoms(self, atoms, translations):
		self.atoms = atoms
		translations = [np.array(t) for t in translations]
		xTot = 0
		yTot = 0
		zTot = 0
		self.xyz = [[],[],[]]
		nat = len(self.atoms)
		self.atomColors = {}
		for atom in self.atoms:
			if atom[0] not in self.atomColors:
				self.atomColors[atom[0]] = self.randomColor()
			vs = geo.onBox(translations, atom[1])
			if len(vs) == 2:
				vs.append(vs[0]+vs[1])
			if len(vs) == 3:
				vs.append(vs[0]+vs[1])
				vs.append(vs[0]+vs[2])
				vs.append(vs[1]+vs[2])
				vs.append(vs[0]+vs[1]+vs[2])
			for t in vs:
				repeat = [x for x in atom]
				repeat[1] = np.array(repeat[1]) + np.array(t)
				self.edgeRepeats.append(repeat)
		for atom in self.atoms+self.edgeRepeats:
			self.xyz[0].append(atom[1][0])
			self.xyz[1].append(atom[1][1])
			self.xyz[2].append(atom[1][2])
			xTot = xTot + atom[1][0]
			yTot = yTot + atom[1][1]
			zTot = zTot + atom[1][2]
		self.cellMaxDim = max([max(self.xyz[0]) - min(self.xyz[0]), max(self.xyz[1]) - min(self.xyz[1]), max(self.xyz[2]) - min(self.xyz[2])])
		self.cellSizeFactor = 0.3/self.cellMaxDim
		self.centerV = 0.5*(translations[0] + translations[1] + translations[2])
		for index, atom in enumerate(self.atoms):
			self.atoms[index][1] = self.atoms[index][1] - self.centerV
		for index, atom in enumerate(self.edgeRepeats):
			self.edgeRepeats[index][1] = self.edgeRepeats[index][1] - self.centerV
		for i, t in enumerate(translations):
			v = t
			M = rz = geo.Rz(-1*geo.xyy(t))
			v = np.inner(rz,v)
			rx = geo.Rx(-1*geo.yzz(v))
			M = np.matmul(rx,M)
			v = np.inner(rx,v)
			ry = geo.Ry(-1*geo.xzz(v))
			M = np.matmul(ry,M)
			self.translations.append([la.norm(t), self.a32l16(la.inv(M)), t])


	def a32l16(self,M):
		return[M[0][0],M[1][0],M[2][0],0,M[0][1],M[1][1],M[2][1],0,M[0][2],M[1][2],M[2][2],0,0,0,0,1]
	def randomColor(self):
				r = random.random();
				g = random.random();
				b = 1 - 0.5*r - 0.5*g
				return (r,g,b,1.)
	def randomizeColors(self):
		for key in self.atomColors:
			self.atomColors[key] = self.randomColor()
		self.updateGL()
	def toggleEdgeAtoms(self, toggle):
		self.xPosDisplay.setText(str(random.random()))


if __name__ == '__main__':
	app = QtGui.QApplication(['Unit Cell'])
	window = MainWindow()
	window.loadAtoms([['Cd',[0.,0.,0.],[1,1,1]],['Cd',[0.5,0.5,0.],[1,1,1]],['Cd',[0.5,0.,0.5],[1,1,1]],['Cd',[0.,0.5,0.5],[1,1,1]],['Te',[0.25,0.25,0.25],[1,1,1]],['Te',[0.75,0.75,0.25],[1,1,1]],['Te',[0.75,0.25,0.75],[1,1,1]],['Te',[0.25,0.75,0.75],[1,1,1]]], [[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
	#window.show()
	window.maximize()
	app.exec_()


