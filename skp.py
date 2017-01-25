import pointgroup as pg
import numpy.linalg as la
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.widgets import Slider, Button, RadioButtons

pi = np.pi

crystal = pg.Crystal()
crystal.addTranslation([1.0,0.0,0.0])
crystal.addTranslation([0.0,1.0,0.0])
crystal.addTranslation([0.0,0.0,1.0])
crystal.addAtom('Cd',[0.0,0.0,0.0])
crystal.addAtom('Cd',[0.5,0.5,0.0])
crystal.addAtom('Cd',[0.5,0.0,0.5])
crystal.addAtom('Cd',[0.0,0.5,0.5])
crystal.addAtom('Te',[0.25,0.25,0.25])
crystal.addAtom('Te',[0.75,0.75,0.25])
crystal.addAtom('Te',[0.75,0.25,0.75])
crystal.addAtom('Te',[0.25,0.75,0.75])
crystal.findPointGroup()
print crystal.system
if not crystal.pointGroup:
	crystal.findPointGroup(checkall=True)
print crystal.pointGroup

kx = 0
ky = 0
kz = 0

M = 10
def filla(M=M, translations=crystal.translations):
	a = []
	v1,v2,v3 = translations
	i = 0
	while len(a) < M:
		if i > 0:
			a.append(i*v1 + i*v2 + i+v3)
		a.append((i+1)*v1 + i*v2 + i*v3)
		a.append(i*v1 + (i+1)*v2 + i*v3)
		a.append(i*v1 + i*v2 + (i+1)*v3)
		a.append((i+1)*v1 + (i+1)*v2 + i*v3)
		a.append((i+1)*v1 + i*v2 + (i+1)*v3)
		a.append(i*v1 + (i+1)*v2 + (i+1)*v3)
		a.append((i+1)*v1 + (i+1)*v2 + (i+1)*v3)
		i += 1
	a = sorted(a, key = lambda v: la.norm(v))
	return a[:M]

a = filla(M, crystal.translations)

def A(a,k,G):
	"""Given a k-vector and a set of symmetry operation matrices, return a plane wave as defined in Evarestov 10.1002/pssb.2221190102 2.14"""
	total = 0
	for g in G:
		total += np.exp(complex(0,1)*np.dot(k,np.matmul(g,a)))
	return total/len(G)


fig = plt.figure()
ax = fig.gca(projection='3d')

# Make data.
X = np.arange(-2*pi, 2*pi, 0.05)
Y = np.arange(-2*pi, 2*pi, 0.05)
X, Y = np.meshgrid(X, Y)
kz = 0.0
surf = None

def makeZ(a):
	Z=[]
	for index, row in enumerate(X):
		Z.append([np.real(A(a,[x,y,kz],crystal.directOps.arrays())) for x,y in zip(row,Y[index])])
	return np.array(Z)

def updateA(val):
	Z = makeZ(np.array([val,val,val]))
	global surf
	surf.remove()
	surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm)
	plt.draw()
	
Z = makeZ([1.0,1.0,0.0])	
#kzAxes = plt.axes([0.25, 0.05, 0.65, 0.03])
#kzSlider = Slider(kzAxes, 'a_m', 0.1, 4.0, valinit=kz, facecolor='k')
#kzSlider.on_changed(updateA)


# Plot the surface.
surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm)

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.show()


















