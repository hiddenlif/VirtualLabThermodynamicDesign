import math,copy,numpy as np
from physics import predict

def number(x,lo,hi):
 if isinstance(x,bool) or not isinstance(x,(int,float)) or not math.isfinite(x) or not lo<=x<=hi:raise ValueError('Numeric action outside domain')
 return float(x)

class Campaign:
 def __init__(self,case):
  self.case=copy.deepcopy(case);self.public=self.case['public'];self.history=copy.deepcopy(case['initial']);self.round=0;self.curves=0;self.stock={'protein_nmol':[10.,10.],'standard_nmol':20.,'ligand_nmol':60.};self.rng=np.random.default_rng(case['noise_seed']);self.done=False
 def state(self):return {'rounds_used':self.round,'curves_used':self.curves,'stock':copy.deepcopy(self.stock)}
 def start(self):return {'type':'start','problem':copy.deepcopy(self.public),'history':copy.deepcopy(self.history),'state':self.state()}
 def batch(self,actions):
  if self.done or self.round>=4 or not isinstance(actions,list) or not 0<=len(actions)<=2:raise ValueError('Batch/round limit')
  stock=copy.deepcopy(self.stock);clean=[]
  for a in actions:
   if not isinstance(a,dict) or set(a)!={'kind','lot','ph','temperature','buffer_h','volume_l','protein_m','syringe_m','injections_l'}:raise ValueError('Action fields')
   if a['kind'] not in ['binding','standard','blank'] or type(a['lot'])is not int or a['lot'] not in [0,1]:raise ValueError('Action kind/lot')
   ph=number(a['ph'],6.5,8.2);T=number(a['temperature'],283.15,310.15);b=number(a['buffer_h'],3600.,47500.)
   if b not in [3600.,21000.,47500.]:raise ValueError('Unknown buffer')
   if number(a['volume_l'],.0002,.0002)!=.0002:raise ValueError('Cell volume')
   p=number(a['protein_m'],5e-6,20e-6);s=number(a['syringe_m'],40e-6,200e-6);inj=a['injections_l']
   if not isinstance(inj,list) or not 4<=len(inj)<=12:raise ValueError('Injection count')
   inject=[number(v,2e-6,15e-6) for v in inj]
   if sum(inject)>.00012+1e-15:raise ValueError('Injection total')
   if a['kind']!='blank':
    stock['protein_nmol'][a['lot']]-=.0002*p*1e9;stock['standard_nmol' if a['kind']=='standard' else 'ligand_nmol']-=sum(inject)*s*1e9
   clean.append(copy.deepcopy(a))
  if min(*stock['protein_nmol'],stock['standard_nmol'],stock['ligand_nmol']) < -1e-9:raise ValueError('Insufficient irreversible stock')
  self.stock=stock;self.round+=1;self.curves+=len(clean);reports=[]
  for a in clean:
   eps=self.public['heat_error'];y=predict(self.case['theta'],a)+self.rng.uniform(-.3,.3,len(a['injections_l']))*eps
   h={'action':a,'y':y.tolist(),'epsilon':eps};reports.append(h);self.history.append(h)
  return {'type':'observations','reports':reports,'state':self.state()}

if __name__ == '__main__':
 import json, sys
 case_file = sys.argv[1] if len(sys.argv) > 1 else 'practice_case.json'
 try:
  with open(case_file) as f:
   case = json.load(f)
  c = Campaign(case)
  msg = c.start()
  print(f"Campaign ready for case with noise seed {case.get('noise_seed')}")
  print(f"Initial stock: protein={c.stock['protein_nmol']}, standard={c.stock['standard_nmol']}, ligand={c.stock['ligand_nmol']}")
 except Exception as e:
  print(f"Error loading campaign: {e}")
  sys.exit(1)
