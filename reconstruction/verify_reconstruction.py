"""Independent numerical cross-checks against preserved historical results, not TAG validation."""
from pathlib import Path
import json,csv,numpy as np
R=Path(__file__).resolve().parent;n=json.loads((R/'results/nominal.json').read_text());s=json.loads((R/'results/sensitivity.json').read_text())
checks=[]
def check(name,ok,details):
 checks.append(dict(check=name,passed=bool(ok),details=details))
 if not ok:raise AssertionError(name)
expected={'triolein_all':608.34,'triolein_measured':653.95,'trilaurin_measured':569.95,'Wu_cap_95':341.1,'Wu_cap_107.9':462.7,'Wu_cap_108.8':441.3,'Wu_cap_109.7':354.4,'Wu_cap_136.9':449.8,'Dymond_full':-103.9,'C17_rotator':260.0,'Wuerflinger_C16_liquid_all':-60.4,'Wuerflinger_C17_liquid_first4':-66.8}
for k,target in expected.items():check('historical_point_'+k,abs(n[k]['E1_MPa']-target)<.055,{'expected_rounded':target,'new':n[k]['E1_MPa'],'tolerance_MPa':.055})
# Independently solve a'(p)+b'(p) dispersion derivative using pairwise lines rather than centred covariance.
for k,v in n.items():
 a=np.array([x['a'] for x in v['isotherms']]);b=np.array([x['b'] for x in v['isotherms']]);da=a[:,None]-a;db=b[:,None]-b;pg=-np.sum(da*db)/np.sum(db*db)
 check('pairwise_identity_'+k,abs(pg-v['E1_MPa'])<1e-8,float(pg))
vs=[x['E1'] for x in s['triolein_measured']['pressure_level_loo']];check('historical_LOO_is_pressure_level',abs(min(vs)-627)<1 and abs(max(vs)-709)<1,[min(vs),max(vs)])
vs2=[x['E1'] for x in s['triolein_measured']['point_loo']];check('point_LOO_different',min(vs2)<min(vs) and max(vs2)>max(vs),[min(vs2),max(vs2)])
t=json.loads((R/'results/temperature_label_sensitivity.json').read_text());check('historical_trilaurin_E2_labels',abs(t['historical_rounded_labels']['E2_MPa']-585.70)<.01,t)
u=json.loads((R/'results/uncertainty_scenarios.json').read_text())
for k,v in u.items():
 for row in v:
  if row['error_model']=='global_scale':check('global_scaling_invariance_'+k+'_'+str(row['relative_sigma']),np.ptp(row['p_quantiles'])<1e-6,row['p_quantiles'])
rr=list(csv.DictReader((R/'data/source_rows.csv').open()));check('no_duplicate_source_rows',len({(x['source'],x['table'],x['substance'],x['group_K'],x['p_MPa']) for x in rr})==len(rr),len(rr))
check('rotator_family_sample_counts',[x['n'] for x in n['C17_rotator']['isotherms']]==[4,4,2],[x['n'] for x in n['C17_rotator']['isotherms']])
j=json.loads((R/'results/C17_junction_volumes.json').read_text());check('source_polynomial_junction',abs(j['p']-239.14)<.01 and abs(j['T']-340.08)<.01,j['p'])
check('volume_jump_reconstruction',abs(j['volumes']['R_to_L']['value_at_junction']-10.75)<.01 and abs(j['volumes']['Cr_to_R']['value_at_junction']-8.01)<.01,j['volumes'])
(R/'results/verification.json').write_text(json.dumps(checks,indent=2));print(len(checks),'checks passed')
