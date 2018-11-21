import numpy as np
import ase.io
from ase.data import covalent_radii
from numpy import linalg as LA

import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D

class Ball:
    def __init__(self,O,R):
        self.O=O
        self.R=R

    def Include(self,P):
        return LA.norm(P-self.O)<=self.R

class PointSet:

    def __init__(self,atoms):
        self.PointSet=[]
        self.maxr=0
        self.n=len(atoms)
        ran=np.random.permutation(self.n)
        for i in range(self.n):
            self.PointSet.append(atoms[ran[i]].position)
            self.maxr=max(self.maxr,covalent_radii[atoms[ran[i]].number])

    def __MakeCircle(self,i,j,k):
        AB=self.PointSet[j]-self.PointSet[i]
        AC=self.PointSet[k]-self.PointSet[i]
        M=np.array([[AB.dot(AB),AB.dot(AC)],[AB.dot(AC),AC.dot(AC)]])
        b=np.array([0.5*AB.dot(AB),0.5*AC.dot(AC)])
        x=LA.solve(M,b)
        AO=x[0]*AB+x[1]*AC
        sbO=self.PointSet[i]+AO
        sbR=LA.norm(sbO-self.PointSet[i])
        return Ball(sbO,sbR)
    
    def __MakeBall(self,i,j,k,l):
        A=self.PointSet[i]
        B=self.PointSet[j]
        C=self.PointSet[k]
        D=self.PointSet[l]
        M=np.array([2*(B-A),2*(C-A),2*(D-A)])
        b=np.array([LA.norm(B)-LA.norm(A),LA.norm(C)-LA.norm(A),LA.norm(D)-LA.norm(A)])
        sbO=LA.solve(M,b)
        sbR=LA.norm(sbO-A)
        return Ball(sbO,sbR) 
        
    def __GetSmallestBallWithijk(self,i,j,k):
        sb=self.__MakeCircle(i,j,k)
        for l in range(k):
            if not sb.Include(self.PointSet[l]):
                tsb=self.__MakeBall(i,j,k,l)
                if tsb.R>sb.R:
                    sb=tsb
        return sb
        
    def __GetSmallestBallWithij(self,i,j):
        sbO=(self.PointSet[i]+self.PointSet[j])/2
        sbR=LA.norm(sbO-self.PointSet[i])
        sb=Ball(sbO,sbR)
        for k in range(j):
            if not sb.Include(self.PointSet[k]):
                sb=self.__GetSmallestBallWithijk(i,j,k)
        return sb
        
    def __GetSmallestBallWithi(self,i):
        sb=Ball(self.PointSet[i],0)
        for j in range(i):
            if not sb.Include(self.PointSet[j]):
                sb=self.__GetSmallestBallWithij(i,j)
        return sb
        
    def GetSmallestBall(self):
        sb=Ball(self.PointSet[0],0)
        for i in range(self.n):
            if not sb.Include(self.PointSet[i]):
                sb=self.__GetSmallestBallWithi(i)
        sb=Ball(sb.O,sb.R+self.maxr)
        return sb

def FindSB(filename):
    atoms=ase.io.read(filename,format='xyz')
    #atoms=ase.io.read(filename,format='vasp')
    a=PointSet(atoms)
    sb=a.GetSmallestBall()
    return sb
    
    
"""
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_aspect("equal")

#atoms=ase.io.read('input',format='xyz')
atoms=ase.io.read('POSCAR',format='vasp')
a=PointSet(atoms)
sb=a.GetSmallestBall()

x=[];y=[];z=[]
for atom in atoms:
    x.append(atom.position[0])
    y.append(atom.position[1])
    z.append(atom.position[2])
ax.scatter(x, y, z,s=a.maxr*100, depthshade=True)

r=sb.R
u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
x=r*np.cos(u)*np.sin(v)+sb.O[0]
y=r*np.sin(u)*np.sin(v)+sb.O[1]
z=r*np.cos(v)+sb.O[2]
ax.plot_wireframe(x, y, z, color="r")

plt.show()
"""
