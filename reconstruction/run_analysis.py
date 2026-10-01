"""Reconstruction dated 2026-09-12. Run python run_analysis.py.
Reads preserved source_rows.csv; no PDF or network required. Not the original August code.
"""
from pathlib import Path
import json,csv,itertools
import numpy as np
from scipy.optimize import minimize_scalar
R=Path(__file__).resolve().parent
(R/'results').mkdir(exist_ok=True)
rows=list(csv.DictReader((R/'data/source_rows.csv').open()))
for x in rows:
 for key in ['T_K','group_K','p_MPa','rho_kg_m3']:
  x[key]=float(x[key])
 x['v_cm3_g']=float(x['v_cm3_g']) if x['v_cm3_g'] else None
def family(pred):
 rr=[x for x in rows if pred(x)];groups={}
 for x in rr:groups.setdefault(x['group_K'],[]).append(x)
 return {t:sorted(xs,key=lambda x:x['p_MPa']) for t,xs in groups.items()}
def fit(f):
 ts=sorted(f);aa=[];bb=[];diag=[]
 for t in ts:
  p=np.array([x['p_MPa'] for x in f[t]]);d=np.array([x['rho_kg_m3'] for x in f[t]])
  b,a=np.polyfit(p,d,1);aa.append(a);bb.append(b);res=d-(a+b*p)
  quad=np.polyval(np.polyfit(p,d,2),p) if len(p)>=3 else np.full_like(p,np.nan)
  diag.append(dict(T=t,n=len(p),p_min=float(p.min()),p_max=float(p.max()),a=a,b=b,rms_linear=float(np.sqrt(np.mean(res**2))),rms_quadratic=(float(np.sqrt(np.mean((d-quad)**2))) if len(p)>=3 else None),residual_signs=''.join('+' if z>0 else '-' for z in res),actual_T_range=[min(x['T_K'] for x in f[t]),max(x['T_K'] for x in f[t])]))
 a=np.array(aa);b=np.array(bb);ac=a-a.mean();bc=b-b.mean();pg=-np.dot(ac,bc)/np.dot(bc,bc)
 at=np.polyfit(ts,a,1)[0];bt=np.polyfit(ts,b,1)[0];e2=-at/bt
 pair=[dict(T1=ts[i],T2=ts[j],p=float((a[j]-a[i])/(b[i]-b[j]))) for i,j in itertools.combinations(range(len(ts)),2) if b[i]!=b[j]]
 return dict(E1_MPa=float(pg),E2_MPa=float(e2),rho_focus=float(np.mean(a+b*pg)),focus_SD_population=float(np.std(a+b*pg)),reach_to_global_max=float(pg/max(x['p_MPa'] for xx in f.values() for x in xx)),isotherms=diag,pairwise=pair,rows=[x['id'] for xx in f.values() for x in xx])
cases={}
base=lambda x:x['source']=='Acosta1996' and x['phase']=='liquid'
for sub in ['triolein','trilaurin']:
 for mode in ['all','measured']:
  cases[f'{sub}_{mode}']=family(lambda x:base(x) and x['substance']==sub and (mode=='all' or x['p_MPa']>.1))
for cap in [95,107.9,108.8,109.7,136.9,999]:cases[f'Wu_cap_{cap}']=family(lambda x:x['source']=='Wu2011' and x['p_MPa']<=cap)
cases['Dymond_full']=family(lambda x:x['source']=='Dymond1979')
cases['C17_rotator']=family(lambda x:x['source']=='Wuerflinger2000' and x['substance']=='C17' and x['phase']=='rotator')
for sub in ['C16','C17']:
 cases[f'Wuerflinger_{sub}_liquid_all']=family(lambda x:x['source']=='Wuerflinger2000' and x['substance']==sub and x['phase']=='liquid')
 cases[f'Wuerflinger_{sub}_liquid_first4']=dict(list(cases[f'Wuerflinger_{sub}_liquid_all'].items())[:4])
results={k:fit(f) for k,f in cases.items()}
# E2 changes with representative temperature labels; E1 does not use T numerically.
temp_variants={}
f=cases['trilaurin_measured'];a=np.array([x['a'] for x in results['trilaurin_measured']['isotherms']]);b=np.array([x['b'] for x in results['trilaurin_measured']['isotherms']])
for name,ts in {'nominal_grid':[328,333,338,343,348,353],'actual_means':[np.mean([x['T_K'] for x in f[t]]) for t in sorted(f)],'first_measured':[f[t][0]['T_K'] for t in sorted(f)],'historical_rounded_labels':[327.6,333,338.1,343.2,348.1,353]}.items():
 temp_variants[name]={'T_labels':ts,'E2_MPa':float(-np.polyfit(ts,a,1)[0]/np.polyfit(ts,b,1)[0])}
(R/'results/temperature_label_sensitivity.json').write_text(json.dumps(temp_variants,indent=2))
# Sensitivity: no optimisation for agreement with historical results.
sens={}
for key in ['triolein_measured','trilaurin_measured','Wu_cap_107.9','Wu_cap_109.7','Dymond_full','C17_rotator']:
 f=cases[key];loo=[]
 for t,xs in f.items():
  if len(xs)<=2:continue
  for i,x in enumerate(xs):
   ff={u:(xx[:i]+xx[i+1:] if u==t else xx) for u,xx in f.items()};loo.append(dict(removed=x['id'],E1=fit(ff)['E1_MPa']))
 famloo=[dict(removed_T=t,E1=fit({u:xx for u,xx in f.items() if u!=t})['E1_MPa']) for t in f if len(f)>2]
 sens[key]=dict(point_loo=loo,isotherm_loo=famloo)
 if key.startswith(('triolein','trilaurin')):
  windows=[]
  for low in [9,24,39,54]:
   for high in [99,114,129,144]:
    ff={t:[x for x in xs if low<=x['p_MPa']<=high] for t,xs in f.items()}
    windows.append(dict(p_min=low,p_max=high,E1=fit(ff)['E1_MPa']))
  sens[key]['windows']=windows
  sens[key]['pressure_level_loo']=[dict(removed_p=p,E1=fit({t:[x for x in xs if x['p_MPa']!=p] for t,xs in f.items()})['E1_MPa']) for p in sorted({x['p_MPa'] for xs in f.values() for x in xs})]
# Parametric uncertainty scenarios: new sensitivity computations, not experimental confidence intervals.
def mc(f,rel,N,seed,mode):
 rng=np.random.default_rng(seed);a=[];b=[]
 global_err=rng.normal(0,rel,(N,1)) if mode=='global_scale' else None
 for t,xs in sorted(f.items()):
  p=np.array([x['p_MPa'] for x in xs]);d=np.array([x['rho_kg_m3'] for x in xs]);X=np.column_stack([np.ones(len(p)),p]);op=np.linalg.pinv(X)
  noise=global_err if mode=='global_scale' else rng.normal(0,rel,(N,1 if mode=='isotherm_scale' else len(p)))
  dd=d*(1+noise);coef=dd@op.T;a.append(coef[:,0]);b.append(coef[:,1])
 a=np.array(a).T;b=np.array(b).T;ac=a-a.mean(1)[:,None];bc=b-b.mean(1)[:,None]
 pg=-np.sum(ac*bc,1)/np.sum(bc**2,1);q=np.quantile(pg,[.025,.5,.975])
 return dict(N=N,seed=seed,relative_sigma=rel,error_model=mode,p_quantiles=q.tolist(),fraction_negative=float(np.mean(pg<0)),note='Scenario interval only; no source-supported independence or complete covariance asserted. Pressure and temperature errors not propagated.')
uncert={}
for i,key in enumerate(['triolein_measured','trilaurin_measured','Wu_cap_107.9','Wu_cap_109.7','Dymond_full','C17_rotator']):
 uncert[key]=[mc(cases[key],rel,20000,20260912+i*100+j*10+m,mode) for j,rel in enumerate([.001,.0025]) for m,mode in enumerate(['point_independent','isotherm_scale','global_scale'])]
# Source representation checks, bounded computational domains explicitly declared.
pgrid=np.linspace(.1,3000,30001);tait={}
for sub,ts,B in [('triolein',[303,313,333,353],[112.73,108.75,100.11,90.76]),('trilaurin',[333,343,353],[91.24,87.69,83.92])]:
 v0=np.array([next(x['v_cm3_g'] for x in rows if x['source']=='Acosta1996' and x['substance']==sub and x['group_K']==t and x['p_MPa']==.1) for t in ts]);B=np.array(B)
 dens=1000/(v0[:,None]*(1-.07325*np.log((pgrid[None,:]+B[:,None])/(.1+B[:,None]))));sd=dens.std(0);imin=int(np.argmin(sd))
 errs=[]
 for j,t in enumerate(ts):
  for x in cases[f'{sub}_measured'][t]:
   v=v0[j]*(1-.07325*np.log((x['p_MPa']+B[j])/(.1+B[j])));errs.append((1000/v-x['rho_kg_m3'])/x['rho_kg_m3'])
 tait[sub]=dict(matched_linear_E1=fit({t:cases[f'{sub}_measured'][t] for t in ts})['E1_MPa'],T=ts,B=B.tolist(),v0=v0.tolist(),c=-.07325,p0=.1,grid_domain=[.1,3000],minimum_p=float(pgrid[imin]),minimum_sd=float(sd[imin]),minimum_at_domain_edge=imin in [0,len(pgrid)-1],strictly_decreasing_on_grid=bool(np.all(np.diff(sd)<0)),fit_MAPD_percent=float(np.mean(np.abs(errs))*100),note='Trilaurin only three Table 19 temperatures; no interpolation of missing B values. Above measured range is mathematical continuation only.')
coef=np.array([[978.3,9.649,0],[926.8,6.773,-.010065],[848.3,5.256,-.003566],[764.8,4.614,-.002066]]);rho0=np.array([770.2,753,735.8,718.6]);p=np.linspace(.1,500,50001);K=coef[:,0,None]+coef[:,1,None]*p+coef[:,2,None]*p*p;d=rho0[:,None]*K/(K-p);sd=d.std(0);j=np.argmin(sd)
eos=dict(p_min=float(p[j]),sd_min=float(sd[j]),domain=[.1,500],coefficients=coef.tolist(),rho0=rho0.tolist(),note='Uses source K=rho*p/(rho-rho0), hence rho=rho0*K/(K-p); extrapolates several branches beyond liquid data. Not terminal evidence.')
# C17 source boundary polynomials and volume jumps are independent calculation tasks.
RL=np.array([294.7,.2347,-1.8784e-4]);CR=np.array([284.6,.2619,-1.2497e-4]);CL=np.array([286.2,.2836,-2.4552e-4]);roots=np.polynomial.polynomial.polyroots(RL-CR);pc=float(roots[roots>0][0]);tc=float(np.polynomial.polynomial.polyval(pc,RL))
volRL=np.array([[.1,30.82],[38.2,27.62],[85.7,23.63],[143.9,18.75]]);volCR=np.array([[.1,9.77],[72.9,9.24],[120.3,8.89],[166,8.55]])
vol={}
for name,arr in [('R_to_L',volRL),('Cr_to_R',volCR)]:
 c=np.polyfit(arr[:,0],arr[:,1],1);vol[name]=dict(input=arr.tolist(),slope=float(c[0]),intercept=float(c[1]),value_at_junction=float(np.polyval(c,pc)),linear_zero=float(-c[1]/c[0]))
junction=dict(R_to_L_coeff=RL.tolist(),Cr_to_R_coeff=CR.tolist(),Cr_to_L_coeff=CL.tolist(),p=pc,T=tc,T_third=float(np.polynomial.polynomial.polyval(pc,CL)),volumes=vol,Cr_to_L_267_6=22.5,excluded_197_7_split=[16.32,9.58],note='Junction recovered from published fitted boundaries, not endpoint-blind prediction. 197.7 MPa partition excluded: source allocation 63:37. Table 4 volumes are reported phase-step determinations, not direct endpoint observations.')
for name,obj in [('nominal',results),('sensitivity',sens),('uncertainty_scenarios',uncert),('source_representations',dict(tait=tait,Dymond= eos)),('C17_junction_volumes',junction)]:
 (R/'results'/f'{name}.json').write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=True))
print('rows',len(rows))
for k,v in results.items():print(k,round(v['E1_MPa'],4),round(v['E2_MPa'],4),round(v['rho_focus'],4),round(v['focus_SD_population'],4))
print('junction',pc,tc,vol);print('representations',tait,eos['p_min'],eos['sd_min'])
