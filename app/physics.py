import math,numpy as np
R=8.31446261815324;T0=298.15;LIGAND=2e-6
THERMAL=(-25000.,-100.,-32000.,-140.)
NAMES=['log_Kint','H_int','Cp_int','pKa_free','pKa_bound','active_0','active_1','heat_offset']
def logk(l,h,c,t):return l+h/R*(1/T0-1/t)+c/R*(np.log(t/T0)+T0/t-1)
def thermo(p,ph,t,b=0):
 p=np.asarray(p);hf,cf,hb,cb=THERMAL
 lf=logk(math.log(10)*p[...,3],hf,cf,t)-math.log(10)*ph;lb=logk(math.log(10)*p[...,4],hb,cb,t)-math.log(10)*ph
 ff=1/(1+np.exp(-lf));fb=1/(1+np.exp(-lb));K=np.exp(logk(p[...,0],p[...,1],p[...,2],t)+np.logaddexp(0,lb)-np.logaddexp(0,lf))
 H=p[...,1]+p[...,2]*(t-T0)+fb*(hb+cb*(t-T0))-ff*(hf+cf*(t-T0))+(fb-ff)*b
 return K,H
def complex(a,b,k):return 2*a*b/(a+b+1/k+np.sqrt((a-b)**2+2*(a+b)/k+1/k**2))
def predict(p,a):
 p=np.asarray(p)
 if a['kind']=='blank':return np.full(len(a['injections_l']),p[7])
 if a['kind']=='standard':k=np.exp(logk(math.log(1e9),-35000.,0.,a['temperature']));h=-35000.
 else:k,h=thermo(p,a['ph'],a['temperature'],a['buffer_h'])
 protein=p[5+a['lot']]*a['protein_m'];ligand=old=0.;out=[]
 for inj in a['injections_l']:
  f=inj/a['volume_l'];protein*=1-f;ligand=ligand*(1-f)+f*a['syringe_m'];new=complex(protein,ligand,k);out.append(a['volume_l']*(new-(1-f)*old)*h+p[7]);old=new
 return np.array(out)
def response(p,lot,ph,t):
 k,_=thermo(p,ph,t);return np.asarray(p)[...,5+lot]*k*LIGAND/(1+k*LIGAND)
def action(kind='binding',lot=0,ph=7.2,t=T0,b=21000.,protein=1e-5,syringe=1e-4):return dict(kind=kind,lot=lot,ph=ph,temperature=t,buffer_h=b,volume_l=.0002,protein_m=protein,syringe_m=syringe,injections_l=[.00001]*8)
