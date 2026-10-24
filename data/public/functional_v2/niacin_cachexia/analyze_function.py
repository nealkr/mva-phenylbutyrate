#!/usr/bin/env python3
"""Source-bound independent-group functional comparisons for two mouse models."""
import argparse,csv,hashlib,itertools,json,math,pathlib,re
import numpy as np
import openpyxl
from scipy import stats

def numeric(value):
    if isinstance(value,(int,float)) and not isinstance(value,bool):
        return float(value),False
    if isinstance(value,str):
        text=value.strip();marked="*" in text
        text=text.replace("*","").replace(",",".")
        if re.fullmatch(r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?",text):
            return float(text),marked
    return None,False

def interval(a,b):
    """Welch 95% CI for treated minus reference mean."""
    a=np.asarray(a,float);b=np.asarray(b,float)
    difference=float(b.mean()-a.mean());va=float(a.var(ddof=1)/len(a));vb=float(b.var(ddof=1)/len(b))
    se=math.sqrt(va+vb)
    if se==0:return difference,[difference,difference],None
    df=(va+vb)**2/(va**2/(len(a)-1)+vb**2/(len(b)-1))
    margin=float(stats.t.ppf(.975,df))*se
    return difference,[difference-margin,difference+margin],float(df)

def permutation(a,b):
    """Two-sided exact mean-difference permutation; preserve endpoint groups."""
    values=np.asarray(list(a)+list(b),float);na=len(a);nb=len(b)
    combinations=math.comb(na+nb,nb)
    if combinations>100000:raise ValueError("Exact permutation exceeds prespecified bound")
    observed=abs(float(np.mean(b)-np.mean(a)));total=float(values.sum());extreme=0
    for selected in itertools.combinations(range(len(values)),nb):
        selected_sum=float(values[list(selected)].sum())
        difference=selected_sum/nb-(total-selected_sum)/na
        extreme+=abs(difference)>=observed-1e-10*max(1,observed)
    return extreme/combinations,combinations

def bh(rows):
    eligible=[x for x in rows if x.get("permutation_p") is not None]
    ordered=sorted(eligible,key=lambda x:x["permutation_p"]);minimum=1.
    for rank in range(len(ordered),0,-1):
        item=ordered[rank-1];minimum=min(minimum,item["permutation_p"]*len(ordered)/rank)
        item["bh_q"]=minimum;item["bh_family_n"]=len(ordered)

def extract(config,book):
    observations=[];omitted=[]
    for endpoint in config["endpoints"]:
        sheet=book[endpoint["sheet"]]
        for group,binding in endpoint["groups"].items():
            cells=[c for row in sheet[binding["range"]] for c in row]
            identities=[c.value for row in sheet[binding["identity_range"]] for c in row] if binding.get("identity_range") else [None]*len(cells)
            baseline=[c for row in sheet[binding["baseline_range"]] for c in row] if binding.get("baseline_range") else [None]*len(cells)
            if len(identities)!=len(cells) or len(baseline)!=len(cells):raise ValueError("Inconsistent within-endpoint identities")
            for cell,identity,base in zip(cells,identities,baseline):
                value,route=numeric(cell.value)
                if base is not None:
                    base_value,base_route=numeric(base.value)
                    value=None if value is None or base_value is None else value-base_value
                    route=route or base_route
                if value is None:
                    if cell.value is not None:
                        if not isinstance(cell.value,str) or not cell.value.strip().startswith("#"):
                            raise ValueError("Unrecognized cell in source observation range: "+cell.coordinate)
                        omitted.append({"endpoint":endpoint["id"],"group":group,"cell":cell.coordinate,"reason":"Published technical exclusion/non-numeric unavailable observation","raw_value":cell.value})
                    continue
                if not math.isfinite(value):raise ValueError("Nonfinite source observation")
                observations.append({"endpoint":endpoint["id"],"model":endpoint["model"],"group":group,"value":value,"unit":endpoint["unit"],"source_sheet":endpoint["sheet"],"source_cell":cell.coordinate,"source_raw_value":str(cell.value),"author_rout_flag":route,"animal_id":identity,"baseline_cell":None if base is None else base.coordinate,"cross_endpoint_pairing_allowed":False})
    return observations,omitted

def analyze(config,observations):
    results=[]
    for policy in ["all_numeric","author_rout_filtered"]:
        for endpoint in config["endpoints"]:
            selected=[x for x in observations if x["endpoint"]==endpoint["id"] and (policy=="all_numeric" or not x["author_rout_flag"])]
            groups={g:[x["value"] for x in selected if x["group"]==g] for g in endpoint["groups"]}
            for reference,treated,contrast in [("disease","treated","niacin_vs_disease"),("control","disease","disease_vs_healthy")]:
                if reference not in groups or treated not in groups:continue
                a,b=groups[reference],groups[treated]
                item={"endpoint":endpoint["id"],"model":endpoint["model"],"role":endpoint["role"],"policy":policy,"contrast":contrast,"unit":endpoint["unit"],"n_reference":len(a),"n_treated":len(b),"author_rout_numeric_n":sum(x["author_rout_flag"] for x in observations if x["endpoint"]==endpoint["id"])}
                if min(len(a),len(b))<3:
                    item["status"]="Insufficient biological observations";results.append(item);continue
                difference,ci,df=interval(a,b);p,nperm=permutation(a,b)
                item.update({"status":"analyzed","reference_mean":float(np.mean(a)),"treated_mean":float(np.mean(b)),"reference_sd":float(np.std(a,ddof=1)),"treated_sd":float(np.std(b,ddof=1)),"mean_difference":difference,"welch_ci95":ci,"welch_df":df,"permutation_p":p,"exact_labelings":nperm})
                if contrast=="niacin_vs_disease" and "control" in groups:
                    healthy=groups["control"]
                    gap=float(np.mean(healthy)-np.mean(a))
                    item["healthy_reference_n"]=len(healthy);item["healthy_mean"]=float(np.mean(healthy))
                    item["fraction_of_healthy_gap_closed"]=None if gap==0 else difference/gap
                    item["healthy_gap_ci95"]=interval(a,healthy)[1]
                    item["gap_closure_denominator_uncertain"]=item["healthy_gap_ci95"][0]<=0<=item["healthy_gap_ci95"][1]
                results.append(item)
        for contrast in ["niacin_vs_disease","disease_vs_healthy"]:
            family=[x for x in results if x["policy"]==policy and x["contrast"]==contrast]
            bh(family)
            for item in family:item["bh_family"]=policy+":"+contrast+":all_inspected_muscle_endpoints_both_models"
    return results

def main():
    parser=argparse.ArgumentParser();parser.add_argument("--config",type=pathlib.Path,default=pathlib.Path(__file__).with_name("computational-config.json"));parser.add_argument("--output",type=pathlib.Path);parser.add_argument("--extract-only",action="store_true");args=parser.parse_args()
    config=json.loads(args.config.read_text());folder=args.config.resolve().parent
    source=folder/config["source"]["file"]
    if hashlib.sha256(source.read_bytes()).hexdigest()!=config["source"]["sha256"]:raise ValueError("Source workbook hash mismatch")
    book=openpyxl.load_workbook(source,data_only=True,read_only=True)
    for check in config["source_checks"]:
        if str(book[check["sheet"]][check["cell"]].value)!=check["value"]:raise ValueError("Source label mismatch")
    observations,omitted=extract(config,book);output=args.output or folder/"prepared";output.mkdir(parents=True,exist_ok=True)
    with (output/"observations.tsv").open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(observations[0]),delimiter="\t");writer.writeheader();writer.writerows(observations)
    inventory={"source_sha256":config["source"]["sha256"],"config_sha256":hashlib.sha256(args.config.read_bytes()).hexdigest(),"n_observations":len(observations),"omitted_non_numeric_cells":omitted,"endpoints":[{"id":e["id"],"model":e["model"],"unit":e["unit"],"role":e["role"],"groups":{g:{"n_numeric":sum(x["endpoint"]==e["id"] and x["group"]==g for x in observations),"n_author_rout":sum(x["endpoint"]==e["id"] and x["group"]==g and x["author_rout_flag"] for x in observations)} for g in e["groups"]}} for e in config["endpoints"]]}
    (output/"source-inventory.json").write_text(json.dumps(inventory,indent=2)+"\n")
    if not args.extract_only:
        result={"schema_version":1,"source_sha256":config["source"]["sha256"],"config_sha256":inventory["config_sha256"],"comparisons":analyze(config,observations),"cohort_reuse":config["cohorts"],"animal_pairing":"Only force's explicit Mouse IDs; no cross-endpoint pairing or correlations","author_known_findings":config["author_known_findings"]}
        (output/"results.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    for path in output.iterdir():
        if path.is_file() and path.stat().st_size>512*1024:raise ValueError("Output table exceeds 512 KiB")
    print(json.dumps({"output":str(output),"endpoints":len(config["endpoints"]),"observations":len(observations),"extract_only":args.extract_only,"files":[{"file":p.name,"bytes":p.stat().st_size,"sha256":hashlib.sha256(p.read_bytes()).hexdigest()} for p in output.iterdir() if p.is_file()]}))

if __name__=="__main__":main()
