"""Independent high-precision replay of complete split-tree certificates."""
import math,time
import enclosures as E
import gradients as G
import summary_bounds as S

def finite(x):
 if isinstance(x,bool) or not isinstance(x,(float,int)) or not math.isfinite(x):raise ValueError('Finite number required')
 return x

def fields(x,names):
 if not isinstance(x,dict) or set(x)!=set(names):raise ValueError('Wrong fields')

def split(box,axis,value):
 if type(axis)is not int or not 0<=axis<len(box):raise ValueError('Invalid split axis')
 value=finite(value)
 if not box[axis][0]<value<box[axis][1]:raise ValueError('Split is not interior')
 a=[b[:] for b in box];b=[b[:] for b in box];a[axis][1]=value;b[axis][0]=value;return a,b

def observation_band(h,k):
 s=E.interval(0);w=E.interval(0)
 for y,inj in zip(h['y'][:k+1],h['action']['injections_l'][:k+1]):
  f=E.interval(inj)/E.interval(h['action']['volume_l']);s=s*(1-f)+E.interval(y);w=w*(1-f)+1
 return s+E.interval(-h['epsilon'],h['epsilon'])*w

def verify(public,history,answer):
 start=time.time();fields(answer,['witness','brackets','determinations','tree']);domain=public['domain'];operating=public['operating'];floor=public['claim_floor']
 witness=answer['witness']
 if not isinstance(witness,list) or len(witness)!=8:raise ValueError('Eight-coordinate witness required')
 for x,b in zip(witness,domain):
  if not b[0]<=finite(x)<=b[1]:raise ValueError('Witness outside prior')
 box=[[x,x] for x in witness]
 for h in history:
  previous=E.interval(0)
  for value,y,inj in zip(E.predict_box(box,h['action']),h['y'],h['action']['injections_l']):
   raw=value-(1-E.interval(inj)/E.interval(h['action']['volume_l']))*previous;previous=value
   band=E.interval(y)+E.interval(-h['epsilon'],h['epsilon'])
   if not (raw.a>=band.a and raw.b<=band.b):raise ValueError('Witness misses observation band')
 brackets=answer['brackets'];labels=answer['determinations']
 if not isinstance(brackets,list) or len(brackets)!=2 or not isinstance(labels,list) or len(labels)!=2:raise ValueError('Two lot results required')
 for (lo,hi),label in zip(brackets,labels):
  lo=finite(lo);hi=finite(hi)
  if not 0<=lo<=hi<=1 or hi-lo>public['max_bracket_width']:raise ValueError('Invalid response bracket')
  if label=='CERTIFIED_COMPLIANT':
   if lo<floor:raise ValueError('Unresolved compliant claim')
  elif label=='CERTIFIED_VIOLATING':
   if hi>=floor:raise ValueError('Unresolved violating claim')
  else:raise ValueError('Unresolved determination')
 counts={'nodes':0,'excluded':0,'retained':0,'operating_nodes':0};witness_retained=False
 def tick(depth,op=False):
  counts['nodes']+=1;counts['operating_nodes']+=int(op)
  if counts['nodes']>public['max_certificate_nodes'] or depth>public['max_certificate_depth']:raise ValueError('Certificate size bound')
 def support_args(proof):
  fields(proof,['channels','weights']);channels=proof['channels'];weights=proof['weights']
  if not isinstance(channels,list) or not 1<=len(channels)<=24 or not isinstance(weights,list) or len(weights)!=len(channels):raise ValueError('Support dimensions')
  for pair in channels:
   if not isinstance(pair,list) or len(pair)!=2:raise ValueError('Support channel')
   j,k=pair
   if type(j)is not int or not 0<=j<len(history) or type(k)is not int or not 0<=k<len(history[j]['y']):raise ValueError('Support channel')
  for w in weights:
   if abs(finite(w))>1e20:raise ValueError('Support weight')
  return channels,weights
 def operating_visit(box,op,node,lot,depth):
  tick(depth,True)
  if node.get('kind')=='monotone':
   fields(node,['kind','support'])
   if not E.decreasing_response(box,op):raise ValueError('Monotonicity not established')
   corner=[[b[1],b[1]] for b in op]
   if node['support'] is None:lo=E.endpoints(E.response_box(box,corner,lot))[0]
   else:
    ch,w=support_args(node['support']);lo=G.support(box,history,ch,w,corner,lot)[0]
   if lo<brackets[lot][0]:raise ValueError('Monotone lower bound fails')
  elif node.get('kind')=='bound':
   fields(node,['kind'])
   if E.endpoints(E.response_box(box,op,lot))[0]<brackets[lot][0]:raise ValueError('Universal lower bound not established')
  elif node.get('kind')=='support':
   fields(node,['kind','channels','weights']);ch,w=support_args({k:node[k] for k in ['channels','weights']})
   if G.support(box,history,ch,w,op,lot)[0]<brackets[lot][0]:raise ValueError('False observation-supported lower bound')
  elif node.get('kind')=='split':
   fields(node,['kind','axis','value','left','right']);a,b=split(op,node['axis'],node['value']);operating_visit(box,a,node['left'],lot,depth+1);operating_visit(box,b,node['right'],lot,depth+1)
  else:raise ValueError('Incomplete operating coverage')
 def visit(box,node,depth):
  nonlocal witness_retained,history
  tick(depth)
  if not isinstance(node,dict):raise ValueError('Malformed tree')
  if node.get('kind')=='calibrate':
   fields(node,['kind','box','blank_ids','standard_ids','child'])
   for key in ['blank_ids','standard_ids']:
    if not isinstance(node[key],list) or any(type(j)is not int or not 0<=j<len(history) for j in node[key]):raise ValueError('Calibration indices')
   new=node['box']
   if not isinstance(new,list) or len(new)!=8:raise ValueError('Calibration dimensions')
   for b in new:
    if not isinstance(b,list) or len(b)!=2:raise ValueError('Calibration coordinate')
    for x in b:finite(x)
   S.validate_calibration(box,history,new,node['blank_ids'],node['standard_ids'])
   visit(new,node['child'],depth+1);return
  if node.get('kind')=='derive':
   fields(node,['kind','constraints','child']);ds=node['constraints']
   if not isinstance(ds,list) or not 1<=len(ds)<=16:raise ValueError('Derived constraint count')
   extra=[]
   for d in ds:
    fields(d,['observation','logK']);j=d['observation'];r=d['logK']
    if type(j)is not int or not 0<=j<len(history):raise ValueError('Derived source')
    if not isinstance(r,list) or len(r)!=2 or not -50<=finite(r[0])<finite(r[1])<=50:raise ValueError('Derived range')
    extra.append(S.derive(box,history,j,r))
   old=history;history=history+extra
   try:visit(box,node['child'],depth+1)
   finally:history=old
   return
  if node.get('kind') in ['linear_contract','linear_exclude']:
   iscontract=node['kind']=='linear_contract';fields(node,['kind','channels','lower_matrix','upper_matrix','box','child'] if iscontract else ['kind','channels','lower_matrix','upper_matrix'])
   ch=node['channels']
   if not isinstance(ch,list) or not 1<=len(ch)<=128:raise ValueError('Linear contraction channels')
   for pair in ch:
    if not isinstance(pair,list) or len(pair)!=2:raise ValueError('Linear contraction channel')
    j,k=pair
    if type(j)is not int or not 0<=j<len(history) or type(k)is not int or not 0<=k<len(history[j]['y']):raise ValueError('Linear contraction index')
   ims=[]
   for key in ['lower_matrix','upper_matrix']:
    mat=node[key]
    if not isinstance(mat,list) or len(mat)!=8:raise ValueError('Linear matrix')
    for row in mat:
     if not isinstance(row,list) or len(row)!=len(ch) or any(abs(finite(x))>1e20 for x in row):raise ValueError('Linear coefficient')
    ims.append(G.image(box,history,mat,ch))
   cut=[[max(a,l[0]),min(b,u[1])] for (a,b),l,u in zip(box,*ims)];empty=any(a>b for a,b in cut)
   if not iscontract:
    if not empty:raise ValueError('False linear exclusion')
    counts['excluded']+=1;return
   new=node['box']
   if empty or not isinstance(new,list) or len(new)!=8:raise ValueError('Linear box')
   for old,claimed,true in zip(box,new,cut):
    if not isinstance(claimed,list) or len(claimed)!=2:raise ValueError('Linear coordinate')
    a,b=map(finite,claimed)
    if not old[0]<=a<=true[0]<=true[1]<=b<=old[1] or not a<b:raise ValueError('Linear contraction omits region')
   visit(new,node['child'],depth+1);return
  if node.get('kind') in ['contract','newton_exclude']:
   iscontract=node['kind']=='contract';fields(node,['kind','channels','matrix','box','child'] if iscontract else ['kind','channels','matrix'])
   channels=node['channels'];matrix=node['matrix']
   if not isinstance(channels,list) or not 1<=len(channels)<=24:raise ValueError('Invalid contraction channel count')
   for pair in channels:
    if not isinstance(pair,list) or len(pair)!=2:raise ValueError('Invalid contraction channel')
    j,k=pair
    if type(j)is not int or not 0<=j<len(history) or type(k)is not int or not 0<=k<len(history[j]['y']):raise ValueError('Invalid contraction channel')
   if not isinstance(matrix,list) or len(matrix)!=8:raise ValueError('Invalid contraction matrix')
   for row in matrix:
    if not isinstance(row,list) or len(row)!=len(channels) or any(abs(finite(c))>1e20 for c in row):raise ValueError('Invalid contraction coefficient')
   im=G.image(box,history,matrix,channels);cut=[[max(a,x),min(b,y)] for (a,b),(x,y) in zip(box,im)]
   empty=any(a>b for a,b in cut)
   if not iscontract:
    if not empty:raise ValueError('False interval-Newton exclusion')
    counts['excluded']+=1;return
   new=node['box']
   if empty or not isinstance(new,list) or len(new)!=8:raise ValueError('Invalid contracted box')
   for old,claimed,true in zip(box,new,cut):
    if not isinstance(claimed,list) or len(claimed)!=2:raise ValueError('Invalid contracted coordinate')
    a,b=map(finite,claimed)
    if not old[0]<=a<=true[0]<=true[1]<=b<=old[1] or not a<b:raise ValueError('Contraction omits compatible region')
   if not any(n!=o for n,o in zip(new,box)):raise ValueError('Empty contraction step')
   visit(new,node['child'],depth+1);return
  if node.get('kind')=='split':
   fields(node,['kind','axis','value','left','right']);a,b=split(box,node['axis'],node['value']);visit(a,node['left'],depth+1);visit(b,node['right'],depth+1)
  elif node.get('kind')=='exclude':
   fields(node,['kind','space','observation','channel']);j=node['observation'];k=node['channel']
   if type(j)is not int or not 0<=j<len(history) or type(k)is not int or not 0<=k<len(history[j]['y']):raise ValueError('Invalid exclusion index')
   h=history[j]
   if node['space']=='cumulative':pred=E.predict_box(box,h['action'])[k];band=observation_band(h,k)
   elif node['space']=='raw':pred=E.predict_raw_box(box,h['action'])[k];band=E.interval(h['y'][k])+E.interval(-h['epsilon'],h['epsilon'])
   else:raise ValueError('Exclusion space')
   if not (pred.b<band.a or pred.a>band.b):raise ValueError('False whole-cell exclusion')
   counts['excluded']+=1
  elif node.get('kind')=='retain':
   fields(node,['kind','lots']);lots=node['lots']
   if not isinstance(lots,list) or len(lots)!=2:raise ValueError('Two retained-lot proofs required')
   for lot,proof in enumerate(lots):
    fields(proof,['lower','upper_point','upper_support']);point=proof['upper_point']
    if not isinstance(point,list) or len(point)!=2 or not all(b[0]<=finite(x)<=b[1] for x,b in zip(point,operating)):raise ValueError('Invalid operating witness')
    if proof['upper_support'] is None:upper=E.endpoints(E.response_box(box,[[x,x] for x in point],lot))[1]
    else:
     ch,w=support_args(proof['upper_support']);upper=G.support(box,history,ch,w,[[x,x] for x in point],lot)[1]
    if upper>brackets[lot][1]:raise ValueError('Worst response upper not established')
    operating_visit(box,operating,proof['lower'],lot,depth+1)
   counts['retained']+=1
   if all(b[0]<=x<=b[1] for x,b in zip(witness,box)):witness_retained=True
  else:raise ValueError('Incomplete parameter coverage')
 visit(domain,answer['tree'],0)
 if not witness_retained:raise ValueError('Witness not in retained region')
 return {'accepted':True,**counts,'seconds':time.time()-start}
