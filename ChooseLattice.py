import numpy as np

def ChooseTriclinic(mins,maxes):
    latticeparm=np.zeros([6])
    latticeparm[0]=np.random.rand()*(maxes[0]-mins[0])+mins[0]
    latticeparm[1]=np.random.rand()*(maxes[1]-mins[1])+mins[1]
    latticeparm[2]=np.random.rand()*(maxes[2]-mins[2])+mins[2]
    latticeparm[3]=np.random.rand()*(maxes[3]-mins[3])+mins[3]
    latticeparm[4]=np.random.rand()*(maxes[4]-mins[4])+mins[4]
    latticeparm[5]=np.random.rand()*(maxes[5]-mins[5])+mins[5]
    return latticeparm

def ChooseMonoclinic(mins,maxes):
    latticeparm=np.zeros([6])
    latticeparm[0]=np.random.rand()*(maxes[0]-mins[0])+mins[0]
    latticeparm[1]=np.random.rand()*(maxes[1]-mins[1])+mins[1]
    latticeparm[2]=np.random.rand()*(maxes[2]-mins[2])+mins[2]
    latticeparm[3]=np.pi/2
    latticeparm[4]=np.random.rand()*(maxes[4]-mins[4])+mins[4]
    latticeparm[5]=np.pi/2
    return latticeparm

def ChooseOrthorhombic(mins,maxes):
    latticeparm=np.zeros([6])
    latticeparm[0]=np.random.rand()*(maxes[0]-mins[0])+mins[0]
    latticeparm[1]=np.random.rand()*(maxes[1]-mins[1])+mins[1]
    latticeparm[2]=np.random.rand()*(maxes[2]-mins[2])+mins[2]
    latticeparm[3]=np.pi/2
    latticeparm[4]=np.pi/2
    latticeparm[5]=np.pi/2
    return latticeparm
    
def ChooseTetragonal(mins,maxes):
    latticeparm=np.zeros([6])
    t1=max(mins[0],mins[1])
    t2=min(maxes[0],maxes[1])
    latticeparm[0]=np.random.rand()*(t2-t1)+t1
    latticeparm[1]=latticeparm[0]
    latticeparm[2]=np.random.rand()*(maxes[2]-mins[2])+mins[2]
    latticeparm[3]=np.pi/2
    latticeparm[4]=np.pi/2
    latticeparm[5]=np.pi/2
    return latticeparm

def ChooseTrigonal(mins,maxes):
    latticeparm=np.zeros([6])
    t1=max(mins[0],mins[1])
    t2=min(maxes[0],maxes[1])
    latticeparm[0]=np.random.rand()*(t2-t1)+t1
    latticeparm[1]=latticeparm[0]
    latticeparm[2]=np.random.rand()*(maxes[2]-mins[2])+mins[2]
    latticeparm[3]=np.pi/2
    latticeparm[4]=np.pi/2
    latticeparm[5]=2*np.pi/3
    return latticeparm

def ChooseHexagonal(mins,maxes):
    latticeparm=np.zeros([6])
    t1=max(mins[0],mins[1])
    t2=min(maxes[0],maxes[1])
    latticeparm[0]=np.random.rand()*(t2-t1)+t1
    latticeparm[1]=latticeparm[0]
    latticeparm[2]=np.random.rand()*(maxes[2]-mins[2])+mins[2]
    latticeparm[3]=np.pi/2
    latticeparm[4]=np.pi/2
    latticeparm[5]=2*np.pi/3 
    return latticeparm

def ChooseCubic(mins,maxes):
    latticeparm=np.zeros([6])
    t1=max(mins[0],mins[1],mins[2])
    t2=min(maxes[0],maxes[1],maxes[2])
    latticeparm[0]=np.random.rand()*(t2-t1)+t1
    latticeparm[1]=latticeparm[0]
    latticeparm[2]=latticeparm[0]
    latticeparm[3]=np.pi/2
    latticeparm[4]=np.pi/2
    latticeparm[5]=np.pi/2
    return latticeparm

def GetLatticeParm(spg,latticeMins,latticeMaxes):
    if spg>=1 and spg<=2:
        latticeparm=ChooseTriclinic(latticeMins,latticeMaxes)
    if spg>=3 and spg<=15:
        latticeparm=ChooseMonoclinic(latticeMins,latticeMaxes)
    if spg>=16 and spg<=74:
        latticeparm=ChooseOrthorhombic(latticeMins,latticeMaxes)
    if spg>=75 and spg<=142:
        latticeparm=ChooseTetragonal(latticeMins,latticeMaxes)
    if spg>=143 and spg<=167:
        latticeparm=ChooseTrigonal(latticeMins,latticeMaxes)
    if spg>=168 and spg<=194:
        latticeparm=ChooseHexagonal(latticeMins,latticeMaxes)
    if spg>=195 and spg<=230:
        latticeparm=ChooseCubic(latticeMins,latticeMaxes)
    ax=latticeparm[0]
    bx=latticeparm[1]*np.cos(latticeparm[5])
    by=latticeparm[1]*np.sin(latticeparm[5])
    cx=latticeparm[2]*np.cos(latticeparm[4])
    cy=(latticeparm[2]*latticeparm[1]*np.cos(latticeparm[3])-cx*bx)/by
    cz=np.sqrt(latticeparm[2]**2-cx**2-cy**2)
    M=np.array([[ax,bx,cx],[0,by,cy],[0,0,cz]])
    return M
    
#M=GetLatticeParm(ChooseCubic,[0,0,0,np.pi/3,np.pi/3,np.pi/3],[5,5,5,2*np.pi/3,2*np.pi/3,2*np.pi/3])