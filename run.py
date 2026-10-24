#!/usr/bin/env python3
"""Run public-source analyses and compare their numerical outputs."""
import argparse,csv,hashlib,json,math,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent
VOLATILE={"generated_utc","python","numpy","scipy"}
def same(a,b,path="root"):
 if isinstance(a,bool) or isinstance(b,bool):
  if a!=b:raise ValueError("Boolean mismatch: "+path)
 elif isinstance(a,(int,float)) and isinstance(b,(int,float)):
  if not math.isclose(a,b,rel_tol=1e-10,abs_tol=1e-10):raise ValueError("Numerical mismatch: "+path)
 elif isinstance(a,dict) and isinstance(b,dict):
  aa={k:v for k,v in a.items() if k not in VOLATILE};bb={k:v for k,v in b.items() if k not in VOLATILE}
  if aa.keys()!=bb.keys():raise ValueError("Object schema mismatch: "+path)
  for k in aa:same(aa[k],bb[k],path+"."+k)
 elif isinstance(a,list) and isinstance(b,list):
  if len(a)!=len(b):raise ValueError("List length mismatch: "+path)
  for i,(x,y) in enumerate(zip(a,b)):same(x,y,path+"["+str(i)+"]")
 elif a!=b:raise ValueError("Value mismatch: "+path)
def check_file(actual,expected):
 if expected.suffix==".json":same(json.loads(actual.read_text()),json.loads(expected.read_text()));return
 with actual.open(newline="") as f:a=list(csv.reader(f,delimiter="\t"))
 with expected.open(newline="") as f:b=list(csv.reader(f,delimiter="\t"))
 if len(a)!=len(b):raise ValueError("Table row-count mismatch")
 for i,(ra,rb) in enumerate(zip(a,b)):
  if len(ra)!=len(rb):raise ValueError("Table column-count mismatch")
  for j,(x,y) in enumerate(zip(ra,rb)):
   try:nx,ny=float(x),float(y)
   except ValueError:same(x,y,str(i)+":"+str(j))
   else:
    if math.isnan(nx) and math.isnan(ny):continue
    same(nx,ny,str(i)+":"+str(j))
def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument("--analysis",action="append",default=[]);ap.add_argument("--core",action="store_true");ap.add_argument("--check",action="store_true");a=ap.parse_args()
 specs=json.loads((ROOT/"analyses.json").read_text())["analyses"];known={s["id"] for s in specs}
 if set(a.analysis)-known:ap.error("Unknown analysis ID")
 selected=[s for s in specs if (a.core and s["group"]=="core") or s["id"] in a.analysis]
 if not selected:ap.error("Choose --core or --analysis using analyses.json")
 receipts=[]
 for s in selected:
  for name in s["reference_files"]:(ROOT/"results"/name).parent.mkdir(parents=True,exist_ok=True)
  for x in s["inputs"]:
   p=ROOT/x["path"]
   if not p.is_file():raise ValueError("Missing input; run fetch_inputs.py first: "+x["path"])
   if hashlib.sha256(p.read_bytes()).hexdigest()!=x["sha256"]:raise ValueError("Input hash mismatch: "+x["path"])
  for command in s["commands"]:
   run=subprocess.run([sys.executable]+command,cwd=ROOT,capture_output=True,text=True)
   if run.returncode:print(run.stderr[-4000:],file=sys.stderr);raise RuntimeError("Analysis failed: "+s["id"])
  compared=[]
  if a.check:
   for name in s["reference_files"]:check_file(ROOT/"results"/name,ROOT/"expected"/name);compared.append(name)
  receipt={"analysis":s["id"],"status":"completed","reference_files_compared":compared};receipts.append(receipt);print(json.dumps(receipt),flush=True)
 (ROOT/"results").mkdir(exist_ok=True)
 (ROOT/"results/reproduction-receipt.json").write_text(json.dumps({"analyses":receipts,"numerical_tolerance":{"relative":1e-10,"absolute":1e-10},"ignored_runtime_metadata_fields":sorted(VOLATILE)},indent=2)+"\n")
if __name__=="__main__":main()
