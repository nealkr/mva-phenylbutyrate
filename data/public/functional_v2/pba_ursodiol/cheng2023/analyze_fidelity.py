#!/usr/bin/env python3
"""Cheng source-percentage analysis; validation-only never evaluates study contrasts."""
import argparse, copy, csv, hashlib, io, itertools, json, math, platform, sys, zipfile
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np
import scipy
from scipy import stats

HASHES = {
 'fidelity_observations.csv':'29e9afc3fe55e03d8744c68a2dbea2bb8cbe2ae0a1843942e0b4765228937727',
 '42255_2023_832_MOESM4_ESM.xlsx':'c7e46081ff1eba40d5cf450ea57bf3202565db1f2ad9432b15dc36ab8d69385d',
 '42255_2023_832_MOESM6_ESM.xlsx':'ed98c767488a673a6833da2593e75325302a321fb5c24d94b28f657f29b8a888'}
HEADERS='file,sheet,cell,endpoint,arm,source_row,within_arm_row_index,value_percent,biological_unit,original_experiment_id,numerator,denominator,pairing_established'.split(',')
ENDPOINTS=['misalignment','lagging_chromosomes','micronuclei']
A2=['DMSO','G6PDi','G6PDi+NAC']
A3=['BAG3 WT','BAG3 WT+NAC','BAG3 T285A','BAG3 T285A+NAC']
CONFIG={
 'schema_version':1,'study_doi':'10.1038/s42255-023-00832-9','input_sha256':HASHES,
 'observation_count':69,'contrast_count':18,'unit':'independent-experiment percentage; contrast in percentage points',
 'endpoints':ENDPOINTS,'pairing':'unpaired; no original IDs or denominators available',
 'simple_contrasts':'Welch t 95% CI and exact two-sided mean-difference label permutation p',
 'interaction':'independent four-group linear contrast with Welch-Satterthwaite 95% CI and t p; no permutation',
 'multiplicity':'Benjamini-Hochberg across all 18 contrasts; simple permutation p plus interaction Welch p',
 'inference_status':'retrospective analysis of inspected published source values; primary flags identify requested rescue contrasts, not preregistration',
 'NAC_endpoint_concentration_mM':None,'NAC_endpoint_duration_hours':None,
 'dose_policy':'Exact Fig2h/3n dose and duration unresolved; studywide 1-5 mM range and 24h imaging window not imputed',
 'confidence_level':0.95,'permutation_alternative':'two-sided absolute unstudentized mean difference',
 'permutation_exchangeability':'Simple label tests assume exchangeable independent observations under the equal-distribution null, stronger than equal means with unequal variances.',
 'primary_rescue_definition':['Fig.2h G6PDi+NAC minus G6PDi','Fig.3n BAG3 T285A+NAC minus BAG3 T285A'],
 'do_not_pool':'Do not pool cells, endpoints, or figures as independent replications; no matched-row or binomial analysis.'}

def require(ok,msg):
 if not ok: raise ValueError(msg)
def digest(b): return hashlib.sha256(b).hexdigest()
def write_json(path,obj): path.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def expected_coordinates():
 d={}
 for fig,fn,arms,cols,maxrows in [('Fig.2h','42255_2023_832_MOESM4_ESM.xlsx',A2,['ABC','EFG','IJK'],[7,7,6]),('Fig.3n','42255_2023_832_MOESM6_ESM.xlsx',A3,['ABCD','FGHI','KLMN'],[6,6,6])]:
  for ep,cs,last in zip(ENDPOINTS,cols,maxrows):
   for arm,col in zip(arms,cs):
    for row in range(4,last+1):d[(fn,fig,col+str(row))]=(ep,arm,row,row-3)
 return d

def source_values(base):
 out={}; ns={'m':'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
 for fn,fig,part in [('42255_2023_832_MOESM4_ESM.xlsx','Fig.2h','xl/worksheets/sheet7.xml'),('42255_2023_832_MOESM6_ESM.xlsx','Fig.3n','xl/worksheets/sheet9.xml')]:
  with zipfile.ZipFile(base/fn) as z:
   require(z.testzip() is None,'ZIP CRC failure')
   root=ET.fromstring(z.read(part))
   for c in root.findall('.//m:sheetData/m:row/m:c',ns):
    key=(fn,fig,c.get('r'))
    if key in expected_coordinates():out[key]=float(c.find('m:v',ns).text)
 return out

def check_rows(rows,src):
 expected=expected_coordinates(); require(len(rows)==69,'row count mismatch');seen=set()
 for r in rows:
  require(list(r)==HEADERS,'schema mismatch')
  key=(r['file'],r['sheet'],r['cell']);require(key in expected,'unexpected source coordinate');require(key not in seen,'duplicate source coordinate');seen.add(key)
  ep,arm,row,index=expected[key]
  require((r['endpoint'],r['arm'],r['source_row'],r['within_arm_row_index'])==(ep,arm,str(row),str(index)),'arm, n, row or endpoint mismatch')
  require(r['biological_unit']=='independent_experiment_per_figure_caption','unit mismatch')
  require(all(r[x]=='' for x in ['original_experiment_id','numerator','denominator']),'invented ID or cell counts')
  require(r['pairing_established']=='False','invented pairing')
  value=float(r['value_percent']);require(math.isfinite(value) and 0<=value<=100,'invalid percentage')
  require(value==src[key],'CSV/source-workbook value mismatch')
 require(seen==set(expected),'missing coordinates')

def validate(base,config_path):
 require(json.loads(config_path.read_text(encoding='utf-8-sig'))==CONFIG,'computational config differs from fixed contract')
 for fn,sha in HASHES.items():require(digest((base/fn).read_bytes())==sha,'SHA mismatch: '+fn)
 rows=list(csv.DictReader(io.StringIO((base/'fidelity_observations.csv').read_text(encoding='utf-8-sig'))))
 src=source_values(base);check_rows(rows,src)
 # Negative checks exercise source gate only, never calculate study outcomes.
 for field,value in [('cell','Z999'),('arm','invented'),('original_experiment_id','invented'),('denominator','100'),('pairing_established','True'),('value_percent','101'),('value_percent','5.91')]:
  bad=copy.deepcopy(rows);bad[0][field]=value
  try:check_rows(bad,src)
  except ValueError:pass
  else:raise AssertionError('negative gate failed: '+field)
 bad=copy.deepcopy(rows);bad[-1]=copy.deepcopy(bad[0])
 try:check_rows(bad,src)
 except ValueError:pass
 else:raise AssertionError('duplicate gate failed')
 require(digest((base/'fidelity_observations.csv').read_bytes()+b' ')!=HASHES['fidelity_observations.csv'],'SHA mutation test failed')
 return rows

def welch_linear(groups,weights):
 require(len(groups)==len(weights),'weight count')
 n=np.array([len(x) for x in groups],dtype=float)
 require(bool(np.all(n>=2)),'n too small')
 means=np.array([np.mean(x) for x in groups]);variances=np.array([np.var(x,ddof=1) for x in groups]);w=np.array(weights,dtype=float)
 est=float(w@means);q=w*w*variances/n;var=float(q.sum());require(var>0,'zero contrast variance')
 se=math.sqrt(var);df=float(var*var/np.sum(q*q/(n-1)));crit=float(stats.t.ppf(.975,df));t=est/se
 return dict(estimate_pp=est,se_pp=se,df=df,ci95_low_pp=est-crit*se,ci95_high_pp=est+crit*se,welch_t=t,welch_p=float(2*stats.t.sf(abs(t),df)))

def exact_permutation(a,b):
 # a minus b; all combinatorial label allocations, including observed allocation.
 pool=np.concatenate([a,b]);na=len(a);total=float(pool.sum());obs=abs(float(np.mean(a)-np.mean(b)));hits=0;count=0
 tol=1e-12*max(1.,obs,float(np.max(np.abs(pool))))
 for ix in itertools.combinations(range(len(pool)),na):
  s=float(pool[list(ix)].sum());stat=abs(s/na-(total-s)/len(b));hits+=int(stat>=obs-tol);count+=1
 return dict(permutation_p=hits/count,permutation_extreme_allocations=hits,permutation_total_allocations=count)

def bh(p):
 p=np.asarray(p,dtype=float);order=np.argsort(p,kind='stable');ranked=p[order]*len(p)/np.arange(1,len(p)+1);adj=np.minimum.accumulate(ranked[::-1])[::-1];out=np.empty(len(p));out[order]=np.minimum(adj,1);return out.tolist()

def toy_tests():
 z=welch_linear([np.array([2.,3.]),np.array([0.,1.])],[1,-1])
 require(abs(z['estimate_pp']-2)<1e-12 and abs(z['se_pp']-math.sqrt(.5))<1e-12 and abs(z['df']-2)<1e-12,'toy Welch test')
 require(abs(exact_permutation(np.array([2.,3.]),np.array([0.,1.]))['permutation_p']-1/3)<1e-12,'toy exact permutation')
 z=welch_linear([np.array([4.,5.]),np.array([2.,3.]),np.array([1.,2.]),np.array([0.,1.])],[1,-1,-1,1])
 require(abs(z['estimate_pp']-1)<1e-12 and abs(z['df']-4)<1e-12 and abs(z['se_pp']-1)<1e-12,'toy interaction test')
 require(np.allclose(bh([.01,.04,.03]),[.03,.04,.04]),'toy BH test')

def analyze(rows):
 result=[]
 for ep in ENDPOINTS:
  for fig,specs in [('Fig.2h',[(A2[1],A2[0],'context',False),(A2[2],A2[1],'rescue',True)]),('Fig.3n',[(A3[2],A3[0],'context',False),(A3[1],A3[0],'reference_NAC',False),(A3[3],A3[2],'rescue',True)])]:
   arms={arm:np.array([float(r['value_percent']) for r in rows if r['sheet']==fig and r['endpoint']==ep and r['arm']==arm]) for arm in (A2 if fig=='Fig.2h' else A3)}
   def append(names,weights,family,primary):
    groups=[arms[a] for a in names];z={'figure':fig,'endpoint':ep,'family':family,'primary_rescue':primary,'arms_in_weight_order':names,'weights':weights,'n':[len(a) for a in groups],'means_percent':[float(np.mean(a)) for a in groups],'sd_percent':[float(np.std(a,ddof=1)) for a in groups]};z.update(welch_linear(groups,weights))
    if len(names)==2:z.update(exact_permutation(*groups));z['bh_input_p']=z['permutation_p'];z['p_method']='exact_unpaired_mean_difference_label_permutation'
    else:z['bh_input_p']=z['welch_p'];z['p_method']='four_group_Welch_Satterthwaite_t'
    result.append(z)
   for a,b,family,primary in specs:append([a,b],[1,-1],family,primary)
   if fig=='Fig.3n':append([A3[3],A3[2],A3[1],A3[0]],[1,-1,-1,1],'NAC_interaction_mutant_minus_WT',False)
 require(len(result)==18,'contrast count failed')
 for z,q in zip(result,bh([x['bh_input_p'] for x in result])):z['bh_q_all18']=q
 return result

def main():
 p=argparse.ArgumentParser();p.add_argument('--data-dir',type=Path,default=Path(__file__).resolve().parent);p.add_argument('--config',type=Path);p.add_argument('--validate-only',action='store_true');p.add_argument('--receipt',type=Path);p.add_argument('--output-dir',type=Path);a=p.parse_args();base=a.data_dir.resolve();cfg=a.config or base/'computational-config.json'
 require(a.validate_only or a.output_dir is not None,'analysis requires --output-dir');require(not(a.validate_only and a.output_dir),'validation mode forbids results output')
 rows=validate(base,cfg);toy_tests()
 receipt={'status':'validated','study_outcomes_computed':not a.validate_only,'observation_count_checked':len(rows),'checks':['fixed_input_SHA256','ZIP_CRC','69_exact_source_coordinates','exact_workbook_value_match','schema_arms_n_percent_range','blank_ID_counts_unpaired','negative_gate_mutations','toy_Welch_exact_permutation_interaction_BH'],'input_sha256':HASHES,'script_sha256':digest(Path(__file__).read_bytes()),'config_sha256':digest(cfg.read_bytes()),'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__}
 if a.receipt:write_json(a.receipt,receipt)
 if a.validate_only:print(json.dumps(receipt,allow_nan=False));return
 a.output_dir.mkdir(parents=True,exist_ok=True);results=analyze(rows);write_json(a.output_dir/'fidelity_results.json',{'config':CONFIG,'validation':receipt,'contrasts':results});write_json(a.output_dir/'run_receipt.json',receipt);print(json.dumps({'status':'complete','contrast_count':len(results),'output':str(a.output_dir)}))

if __name__=='__main__':main()