"""Interval automatic differentiation of cumulative thermodynamic heats."""
import enclosures as E
N=8
class D:
 def __init__(self,v,g=None):self.v=E.interval(v) if isinstance(v,(int,float)) else v;self.g=[E.interval(0) for _ in range(N)] if g is None else g
 def __add__(self,b):b=cv(b);return D(self.v+b.v,[a+c for a,c in zip(self.g,b.g)])
 __radd__=__add__
 def __neg__(self):return D(-self.v,[-a for a in self.g])
 def __sub__(self,b):return self+-cv(b)
 def __rsub__(self,b):return cv(b)+-self
 def __mul__(self,b):b=cv(b);return D(self.v*b.v,[a*b.v+self.v*c for a,c in zip(self.g,b.g)])
 __rmul__=__mul__
 def __truediv__(self,b):b=cv(b);return D(self.v/b.v,[(a*b.v-self.v*c)/(b.v*b.v) for a,c in zip(self.g,b.g)])
 def __rtruediv__(self,b):return cv(b)/self
 def __pow__(self,n):
  if n!=2:raise ValueError('power')
  return D(self.v**2,[2*self.v*g for g in self.g])
def cv(x):return x if isinstance(x,D) else D(x)
def exp(x):x=cv(x);v=E.iv.exp(x.v);return D(v,[v*g for g in x.g])
def ln(x):x=cv(x);return D(E.iv.ln(x.v),[g/x.v for g in x.g])
def sigmoid(x):
 x=cv(x);lo,hi=E.endpoints(x.v)
 def val(z):return 1/(1+E.iv.exp(-E.interval(z)))
 a,b=E.endpoints(val(lo))[0],E.endpoints(val(hi))[1];v=E.interval(a,b)
 da=E.interval(a)*(1-E.interval(a));db=E.interval(b)*(1-E.interval(b))
 low=min(E.endpoints(da)[0],E.endpoints(db)[0]);high=.25 if a<=.5<=b else max(E.endpoints(da)[1],E.endpoints(db)[1])
 return D(v,[E.interval(low,high)*g for g in x.g])
def softplus(x):
 x=cv(x);v=E.iv.ln(1+E.iv.exp(x.v));factor=sigmoid(x).v
 return D(v,[factor*g for g in x.g])
def sqrt(x):x=cv(x);v=E.iv.sqrt(x.v);return D(v,[g/(2*v) for g in x.g])
def logk(l,h,c,t):
 # Temperature here is a measured point, so no operating-variable derivative.
 return l+h/E.R*(1/E.T0-1/E.interval(t))+c/E.R*(E.iv.ln(E.interval(t)/E.T0)+E.T0/E.interval(t)-1)
def predict(box,a):
 p=[]
 for i,b in enumerate(box):g=[E.interval(0) for _ in range(N)];g[i]=E.interval(1);p.append(D(E.interval(*b),g))
 if a['kind']=='derived':
  inner=dict(a['source']);inner['kind']='thermo_only';values=predict(box,inner);return [(v-a['centers'][i])/a['scales'][i] for i,v in enumerate(values)]
 if a['kind']=='blank':k=h=None
 elif a['kind']=='standard':k=cv(E.iv.exp(logk(cv(E.iv.ln(E.interval(1e9))),-35000,0,a['temperature']).v));h=cv(-35000)
 else:
  t=a['temperature'];ln10=E.iv.ln(E.interval(10));lf=logk(p[3]*ln10,-25000,-100,t)-ln10*a['ph'];lb=logk(p[4]*ln10,-32000,-140,t)-ln10*a['ph'];ff=sigmoid(lf);fb=sigmoid(lb)
  logK=logk(p[0],p[1],p[2],t)+softplus(lb)-softplus(lf);k=exp(logK);dt=E.interval(t)-E.T0
  h=p[1]+p[2]*dt+fb*(a['buffer_h']-32000-140*dt)-ff*(a['buffer_h']-25000-100*dt)
  if a['kind']=='thermo_only':return [logK,h]
 v=a['volume_l'];protein=p[5+a['lot']]*a['protein_m'];ligand=E.interval(0);gain=E.interval(0);out=[]
 for inj in a['injections_l']:
  f=E.interval(inj)/E.interval(v);gain=gain*(1-f)+1
  if a['kind']=='blank':out.append(p[7]*gain);continue
  protein=protein*(1-f);ligand=ligand*(1-f)+f*a['syringe_m'];b=cv(ligand)
  s=protein+b+1/k;disc=(protein-b)**2+2*(protein+b)/k+1/(k*k);complex=2*protein*b/(s+sqrt(disc));out.append(v*complex*h+p[7]*gain)
 return out

def image(box,history,matrix,channels):
 """Mean-value/Krawczyk enclosure valid for every band-compatible point."""
 mid=[a+(b-a)/2 for a,b in box];center=[[x,x] for x in mid];delta=[E.interval(a,b)-E.interval(x) for (a,b),x in zip(box,mid)];values=[];jac=[];bands=[]
 cached={}
 for obs,k in channels:
  h=history[obs]
  if obs not in cached:cached[obs]=(predict(center,h['action']),predict(box,h['action']))
  at,wide=cached[obs];s=E.interval(0);w=E.interval(0)
  for y,inj in zip(h['y'][:k+1],h['action']['injections_l'][:k+1]):
   f=E.interval(inj)/E.interval(h['action']['volume_l']);s=s*(1-f)+y;w=w*(1-f)+1
  bands.append(s+E.interval(-h['epsilon'],h['epsilon'])*w);values.append(at[k].v);jac.append(wide[k].g)
 out=[]
 for i,row in enumerate(matrix):
  z=E.interval(mid[i])-sum((E.interval(c)*value for c,value in zip(row,values)),E.interval(0))
  # Combine coefficients before bounding shared raw measurement errors.
  coeff={}
  for c,(obs,k) in zip(row,channels):
   weight=E.interval(c);h=history[obs]
   for q in range(k,-1,-1):
    coeff[obs,q]=coeff.get((obs,q),E.interval(0))+weight
    weight*=1-E.interval(h['action']['injections_l'][q])/E.interval(h['action']['volume_l'])
  for (obs,k),c in coeff.items():
   h=history[obs];z+=c*(E.interval(h['y'][k])+E.interval(-h['epsilon'],h['epsilon']))
  for j in range(8):z+=(E.interval(int(i==j))-sum((E.interval(c)*g[j] for c,g in zip(row,jac)),E.interval(0)))*delta[j]
  out.append(E.endpoints(z))
 return out

def response(box,op,lot):
 p=[]
 for i,b in enumerate(box):g=[E.interval(0) for _ in range(N)];g[i]=E.interval(1);p.append(D(E.interval(*b),g))
 ph=E.interval(*op[0]);t=E.interval(*op[1]);aa=(1/E.T0-1/t)/E.R;bb=E.logk(E.interval(0),E.interval(0),E.interval(1),t)
 def lk(l,h,c):return cv(l)+cv(h)*aa+cv(c)*bb
 ln10=E.iv.ln(E.interval(10));lf=lk(p[3]*ln10,-25000,-100)-ln10*ph;lb=lk(p[4]*ln10,-32000,-140)-ln10*ph
 logK=lk(p[0],p[1],p[2])+softplus(lb)-softplus(lf)
 return p[5+lot]*sigmoid(logK+E.iv.ln(E.LIGAND))

def support(box,history,channels,weights,op,lot):
 mid=[a+(b-a)/2 for a,b in box];center=[[x,x] for x in mid];delta=[E.interval(a,b)-E.interval(x) for (a,b),x in zip(box,mid)]
 f0=response(center,op,lot).v;df=response(box,op,lot).g;cached={};values=[];jac=[]
 for j,k in channels:
  if j not in cached:cached[j]=(predict(center,history[j]['action']),predict(box,history[j]['action']))
  a,b=cached[j];values.append(a[k].v);jac.append(b[k].g)
 z=f0-sum((E.interval(w)*v for w,v in zip(weights,values)),E.interval(0));coeff={}
 for c,(j,k) in zip(weights,channels):
  w=E.interval(c);h=history[j]
  for q in range(k,-1,-1):coeff[j,q]=coeff.get((j,q),E.interval(0))+w;w*=1-E.interval(h['action']['injections_l'][q])/E.interval(h['action']['volume_l'])
 for (j,k),w in coeff.items():h=history[j];z+=w*(E.interval(h['y'][k])+E.interval(-h['epsilon'],h['epsilon']))
 for i in range(8):z+=(df[i]-sum((E.interval(w)*row[i] for w,row in zip(weights,jac)),E.interval(0)))*delta[i]
 return E.endpoints(z)

# Reuse arithmetic within one certificate cell; no optimizer or hidden data cache.
from functools import lru_cache
import json
_predict=predict
@lru_cache(maxsize=128)
def _cached_predict(box,encoded):return _predict([list(b) for b in box],json.loads(encoded))
def predict(box,action):return _cached_predict(tuple(tuple(b) for b in box),json.dumps(action,sort_keys=True,separators=(',',':')))
