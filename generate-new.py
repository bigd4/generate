import numpy as np
import copy
from ChooseLattice import GetLatticeParm
import yaml
import os
import ase.io
import sb
from ase.data import covalent_radii, atomic_numbers
from numpy import linalg as LA

def Standford(p):
    p=p-p.astype(int)+1
    return p-p.astype(int)

def CalDistance(p1,p2,latticeparm):
    return np.sum((p1.dot(latticeparm)-p2.dot(latticeparm))**2)

def CheckDistance(p1s,p2s,r1,r2,latticeparm,threshold):
    for i in p1s:
        for j in p2s:
            if CalDistance(i,j,latticeparm)<((r1+r2)*threshold)**2:
                return(False)
    return(True)  

def str_(l):
    l=str(l)
    return l[l.find('[')+1:l.find(']')]+'\n'

class Atom:
    def __init__(self, name, radius=0, module=False, position=None):
        self.name = name
        self.radius=radius
        if module:
            self.PointSet=[]
            b=sb.FindSB(name)
            self.O=b.O
            if radius==0:
                self.radius=b.R

            atoms=ase.io.read(name,format='xyz')
            for atom in atoms:
                self.PointSet.append(Atom(name=atom.symbol,position=atom.position))
        if (not module) and (radius==0):
            self.radius=covalent_radii[atomic_numbers[self.name]]
        if position is not None:
            self.position=position
class Atoms:
    def __init__(self,atom,number,wycks):
        self.atom=atom
        self.number=number
        self.left=number
        self.wycks=wycks
        self.positions=[]
    
    def RemoveDuplication(self,positions,latticeparm,threshold):
        selected=[]
        for position1 in positions:
            shouldadd=True
            for position2 in selected:
                if CalDistance(position1,position2,latticeparm)<(2*self.atom.radius*threshold)**2:
                    shouldadd=False
                    break
            if shouldadd:
                selected.append(position1)    
        return(selected)

    def AddAtom(self,atom):
        self.number+=1
        self.positions.append(atom.position)

class Structure:
    def __init__(self, structure):
        self.structure=structure
        #self.legal=True

    def AllAtomsUsed(self):
        _=np.array([atoms.left for atoms in self.structure])
        return np.max(_>0)==0
    
    def AtomsCanUse(self):
        for label,atoms in enumerate(self.structure):
            if atoms.left>0:
                return label   

    def AddAtoms(self,atoms):
        self.structure.append(atoms)

    def AddWyck(self,atomslabel,wycks,wycklabel):
        structure_=copy.deepcopy(self.structure)
        wycks=copy.deepcopy(wycks)
        structure_[atomslabel].left-=wycks[wycklabel].multiplicity
        if wycks[wycklabel].unique==True:
            wycks[wycklabel].usable=False
        structure_[atomslabel].wycks.append(wycks[wycklabel])
        return Structure(structure_),wycks

    def getVolume(self):
        self.V=np.linalg.det(self.latticeparm)
        return self.V

    def ChooseLattice(self,latticeMins,latticeMaxes,volumeMin,volumeMax):
        volume=0
        while volume>volumeMax or volume<volumeMin:
            self.latticeparm=GetLatticeParm(self.spg,latticeMins,latticeMaxes)
            volume=self.getVolume()

    def MakeCrystal(self):
        num=0
        for i,atoms in enumerate(self.structure):
            for wyck in atoms.wycks:
                attempt=0
                while attempt<wyck.multiplicity*num*5+200:
                    shouldadd=True
                    for j in range(100):
                        positions=atoms.RemoveDuplication(wyck.GetAllPosition(wyck.GetOnePosition()),self.latticeparm,threshold=1)
                        if len(positions)==wyck.multiplicity:
                            break   
                    if len(positions)<wyck.multiplicity:
                        self.legal=False
                        return self
                    for j in range(i):
                        if not CheckDistance(positions,self.structure[j].positions,atoms.atom.radius,self.structure[j].atom.radius,self.latticeparm,threshold=1):
                            shouldadd=False
                            break
                    if shouldadd:
                        self.structure[i].positions.extend(positions)
                        num=num+len(positions)
                        break
                    attempt+=1
                if not shouldadd:
                    self.legal=False
                    return self    
        self.legal=True        
        return self
    
            
    def Displace(self):
        news=Structure([])
        news.legal=True
        news.latticeparm=self.latticeparm
        for atoms in self.structure:
            if 'module' not in atoms.atom.name:
                news.AddAtoms(atoms)
            
        for atoms in self.structure:
            if 'module' in atoms.atom.name:
                for position in atoms.positions:
                    theta=2*np.random.rand()*np.pi
                    phi=2*np.random.rand()*np.pi
                    for atom in atoms.atom.PointSet:
                        Exist=False
                        a=Atom(atom.name,position=self.newposition(position,atom.position,atoms.atom.O,theta,phi))
                        for atoms_ in news.structure:
                            if atoms_.atom.name==atom.name:              
                                atoms_.AddAtom(a)
                                Exist=True
                                break
                        if not Exist:
                            news.AddAtoms(Atoms(Atom(atom.name),0,[]))
                            news.structure[-1].AddAtom(a) 
        return news
    
        
    def newposition(self,p0,p1,O,theta,phi):
        p0=np.array(p0)
        p1=np.array(p1)
        O=np.array(O)
        vec=p1-O
        M1=np.array([[np.cos(theta),-np.sin(theta),0],[np.sin(theta),np.cos(theta),0],[0,0,1]])
        M2=np.array([[1,0,0],[0,np.cos(phi),-np.sin(phi)],[0,np.sin(phi),np.cos(phi)]])
        newvec=M1.dot(M2).dot(vec)
        fvec=LA.solve(self.latticeparm.T,newvec)    
        return p0+fvec
        
    def WritePoscar(self,filename):
        f=open(filename,'w')
        f.write('{}\n'.format(self.legal))
        if self.legal:
            f.write('1 \n')
            f.write(str_(self.latticeparm[0]))
            f.write(str_(self.latticeparm[1]))
            f.write(str_(self.latticeparm[2]))
            for atoms in self.structure:
                f.write(str(atoms.atom.name)+'   ')
            f.write('\n')
            for atoms in self.structure:
                f.write(str(atoms.number)+'   ')
            f.write('\nDirect\n') 
            for atoms in self.structure:
                for position in atoms.positions:
                        f.write(str_(position))
        f.close()


class WyckoffPosition:
    
    def __init__(self,label,wyckmatrix,rotatematrix,transmatrix,multiplicity,unique):
        self.label=label
        self.wyckmatrix=wyckmatrix
        self.rotatematrix=rotatematrix
        self.transmatrix=transmatrix
        self.unique=unique
        self.multiplicity=multiplicity
        self.usable=True

    def IsUsable(self,atoms):
        if self.multiplicity>atoms.left or self.usable==False:
            return False
        if atoms.wycks!=[]:
            if self.label<atoms.wycks[-1].label:
                return False
        return True

    def GetOnePosition(self):
        return Standford(self.wyckmatrix[:,:-1].dot(np.random.rand(3))+self.wyckmatrix[:,-1])
    
    def GetAllPosition(self,p):
        tempositions=[]
        positions=[]
        for rotate in self.rotatematrix:
            p=rotate[:,:-1].dot(p)+rotate[:,-1]
            tempositions.append(p)
    
        for translate in self.transmatrix:
            for position in tempositions:
                positions.append(Standford(translate+position))
        return(positions)  

        

def GetAllCombination(structure,wycks,combinations):
    if structure.AllAtomsUsed():
        combinations.append(copy.deepcopy(structure))
        return combinations
    atomslabel=structure.AtomsCanUse()
    for wycklabel,wyck in enumerate(wycks):
        if wyck.IsUsable(structure.structure[atomslabel]):
            [structure_,wycks_]= structure.AddWyck(atomslabel,wycks,wycklabel)
            structure_.spg=structure.spg
            combinations=GetAllCombination(structure_,wycks_,combinations)
    return combinations

def Initialize(spg,atomlist,installdir):
    wycks=[]
    trans=yaml.load(open(installdir+'/trans-mat/trans'+str(spg)+'.yaml'))
    r=np.array(trans['symmetry'])
    t=np.array(trans['translate'])
    info=yaml.load(open(installdir+'/wyckoff-mat/spg'+str(spg)+'.yaml'))
    for wyck in info:
        l=wyck['label']
        m=np.array(wyck['matrix'])
        n=wyck['multiplicity']
        u=wyck['unique']
        wycks.append(WyckoffPosition(l,m,r,t,n,u))
    structure=[]
    atomvolume=0
    maxr=0
    for atom in atomlist:
        structure.append(Atoms(Atom(atom[0],atom[1],module=atom[2]),atom[3],[]))
        atomvolume+=4*np.pi/3*structure[-1].number*(structure[-1].atom.radius)**3
        maxr=max(maxr,structure[-1].atom.radius)
    s=Structure(structure)
    s.atomvolume=atomvolume
    s.spg=spg
    s.maxr=maxr
    return wycks,s

def CompleteParm(info,s):
    pwd=os.getcwd()
    if 'minVolume' not in info.keys():
        info['minVolume']=s.atomvolume*0.5
    if 'maxVolume' not in info.keys():
        info['maxVolume']=s.atomvolume*2
    if 'outdir' not in info.keys():
        info['anspath']=pwd+'/out/'
    else:
        info['anspath']=pwd+'/{}/'.format(info['outdir'])
    if 'threshold' not in info.keys():
        info['threshold']=1
    if 'latticeMins' not in info.keys():
        if 2*s.maxr<3:
            info['latticeMins']=[3.0,3.0,3.0,60.0,60.0,60.0]
            maxr=3
        else:
            info['latticeMins']=[2*s.maxr,2*s.maxr,2*s.maxr,60.0,60.0,60.0]
            maxr=2*s.maxr
    maxlen=info['maxVolume']/(maxr**2)
    if 'latticeMaxes' not in info.keys():
        if maxlen>3:
            info['latticeMaxes']=[maxlen,maxlen,maxlen,120.0,120.0,120.0]
        else:
            info['latticeMaxes']=[4.0,4.0,4.0,120.0,120.0,120.0]
    return info


def Generate(info):
    atomlist=info['atomlist']
    spg=info['spg']
    installdir=info['installdir']
    [wycks,structure]=Initialize(spg,atomlist,installdir)
    info=CompleteParm(info,structure)
    maxAttempts=info['maxAttempts']
    minVolume=info['minVolume']
    maxVolume=info['maxVolume']
    latticeMins=info['latticeMins']
    latticeMaxes=info['latticeMaxes']
    
    ans=[]
    for i in range(3,6):
        latticeMins[i]=latticeMins[i]*np.pi/180
        latticeMaxes[i]=latticeMaxes[i]*np.pi/180
    
    
    combinations=GetAllCombination(structure,wycks,[])
    if len(combinations)==0:
        return structure,'Fail to find combinations'
    while len(ans)<info['spgnumber']:
        for j_ in range(maxAttempts):
            s=copy.deepcopy(combinations[np.random.randint(0,len(combinations))])
            s.ChooseLattice(latticeMins,latticeMaxes,minVolume,maxVolume)
            s.MakeCrystal()
            if s.legal==True:
                break
        if s.legal==True:
            s=s.Displace()
            ans.append(s)
    return ans,'Success'

d={
    'atomlist':[['Mg',0,False,2],['module_N5',2.5,True,4]],'spg':1,'spgnumber':1,
    'maxAttempts':1000,'installdir':'/fs08/home/js_wangjj/generate-new',
    #'latticeMaxes':[100,100,100, 120.0, 120.0, 120.0],                                                                                                             
    #'latticeMins':[10,10,10, 60.0, 60.0, 60.0],
    #'maxVolume':8003.148547966757,'minVolume':2000.78713699168924
}
[ans,suc]=Generate(d)
if suc=='Success':
    for i,s in enumerate(ans):
        s.WritePoscar('/fs08/home/js_wangjj/generate-new/rua/N5_{}.vasp'.format(i))
