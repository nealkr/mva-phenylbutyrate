#!/usr/bin/env python3
"""Acquire only the version-bound public inputs named in provenance.json."""
import argparse,csv,hashlib,json,shutil,subprocess,sys,urllib.request
from pathlib import Path
ROOT=Path(__file__).resolve().parent
def sha(path):
 h=hashlib.sha256()
 with path.open("rb") as f:
  for b in iter(lambda:f.read(1024*1024),b""):h.update(b)
 return h.hexdigest()
def local(name):
 p=ROOT/name
 if p.is_symlink() or not p.resolve().is_relative_to(ROOT):raise ValueError("Input path escapes release")
 return p
def acquire(s,cache):
 p=local(s["path"])
 if p.exists():
  if sha(p)!=s["sha256"]:raise ValueError("Existing input hash differs: "+s["path"])
  return
 p.parent.mkdir(parents=True,exist_ok=True)
 if cache:
  q=cache/s["path"]
  if q.is_file() and sha(q)==s["sha256"]:shutil.copyfile(q,p);return
 part=p.with_suffix(p.suffix+".part");errors=[]
 for url in s.get("urls",[]):
  try:
   req=urllib.request.Request(url,headers={"User-Agent":"Public-analysis-reproduction/1.0"});count=0
   with urllib.request.urlopen(req,timeout=90) as response,part.open("wb") as f:
    while block:=response.read(65536):
     count+=len(block)
     if count>s["bytes"]:raise ValueError("Response exceeds frozen size")
     f.write(block)
   if count!=s["bytes"] or sha(part)!=s["sha256"]:raise ValueError("Response differs from frozen source")
   part.replace(p);return
  except Exception as e:errors.append(type(e).__name__+": "+str(e));part.unlink(missing_ok=True)
 raise RuntimeError("Could not fetch "+s["id"]+"; manually obtain the exact file in provenance.json. "+"; ".join(errors))
def predictions():
 import pandas as pd
 d=ROOT/"data/public/bubr1-functional-v2"
 with (d/"prediction-labels.tsv").open(newline="") as f:rows=list(csv.DictReader(f,delimiter="\t"))
 source=pd.read_csv(d/"rasp-O60566-v2.csv").set_index("variant")
 if not source.index.is_unique:raise ValueError("Duplicate RaSP identity")
 for r in rows:
  r["rasp_score_ml"]=str(float(source.loc[r["variant"],"score_ml"]))
  r["source_model_plddt"]=str(float(source.loc[r["variant"],"b_factors"]))
 fields=["variant","residue","reduced_abundance_label","expanded_only","measurement_source","alphamissense","rasp_score_ml","source_model_plddt"]
 with (d/"prediction-controls.tsv").open("w",newline="") as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter="\t",lineterminator="\n");w.writeheader();w.writerows(rows)
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument("--group",choices=["core","import","all"],default="core");ap.add_argument("--cache-dir",type=Path);a=ap.parse_args()
 m=json.loads((ROOT/"provenance.json").read_text())
 for s in m["inputs"]:
  if a.group=="all" or a.group in s["groups"]:acquire(s,a.cache_dir);print(json.dumps({"input":s["id"],"sha256_verified":True}))
 if a.group in ["core","all"]:
  predictions()
  subprocess.run([sys.executable,"data/public/functional_v2/pba_ursodiol/cheng2023/extract_fidelity.py"],cwd=ROOT,check=True,stdout=subprocess.DEVNULL)
 for s in m["derived_inputs"]:
  if (a.group=="all" or a.group in s["groups"]) and sha(local(s["path"]))!=s["sha256"]:raise ValueError("Derived input differs: "+s["path"])
 print(json.dumps({"group":a.group,"status":"inputs_verified"}))
if __name__=="__main__":main()
