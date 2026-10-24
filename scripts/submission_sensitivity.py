#!/usr/bin/env python3
"""Fixed public PBA deletion and free-exposure requirement sensitivity."""
import argparse
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path
from statistics import fmean
from openpyxl import load_workbook


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def exact_test(a, b):
    delta, values = fmean(a)-fmean(b), a+b
    extreme = total = 0
    for indices in itertools.combinations(range(len(values)), len(a)):
        chosen = set(indices)
        d = fmean(v for i,v in enumerate(values) if i in chosen)-fmean(v for i,v in enumerate(values) if i not in chosen)
        extreme += abs(d) >= abs(delta)-1e-12
        total += 1
    return delta, extreme/total, extreme, total


def write_tsv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    ap.add_argument("--design", default="results/submission_sensitivity/design.json")
    ap.add_argument("--out", default="results/submission_sensitivity")
    args = ap.parse_args()
    root, design_path = args.root.resolve(), args.root.resolve()/args.design
    design = json.loads(design_path.read_text())
    paths = {k:root/v["path"] for k,v in design["inputs"].items()}
    for k,p in paths.items():
        if digest(p) != design["inputs"][k]["sha256"]:
            raise ValueError(f"Frozen public input changed: {k}")
    cfg = json.loads(paths["neuronal_config"].read_text())
    baseline = json.loads(paths["neuronal_summary"].read_text())
    assert baseline["config_sha256"] == digest(paths["neuronal_config"])
    assert baseline["source_sha256"] == digest(paths["neuronal_workbook"])
    with paths["neuronal_observations"].open() as stream:
        deposited = {(r["outcome"],r["source_cell"]):float(r["value"]) for r in csv.DictReader(stream,delimiter="\t")}
    book = load_workbook(paths["neuronal_workbook"], data_only=True, read_only=True)
    sheet, rows, neuronal = book[cfg["sheet"]], [], []
    assert [r["id"] for r in cfg["outcomes"]] == design["neuronal_methods"]["endpoints"]
    for endpoint in cfg["outcomes"]:
        name, arms = endpoint["id"], {}
        assert sheet[endpoint["panel_cell"]].value == endpoint["panel"]
        for group in ("treated","vehicle"):
            cells = [f"{cfg['groups'][group]['column']}{r}" for r in endpoint["data_rows"]]
            values = [float(sheet[c].value) for c in cells]
            assert len(values)==3 and all(math.isfinite(v) and v>=0 for v in values)
            assert all(math.isclose(v,deposited[name,c],abs_tol=1e-12) for c,v in zip(cells,values))
            arms[group] = list(zip(cells,values))
        subsets = [("full",None,None)]
        subsets += [("delete_one_single_arm",i,None) for i in range(3)]
        subsets += [("delete_one_single_arm",None,j) for j in range(3)]
        subsets += [("delete_one_each_arm",i,j) for i in range(3) for j in range(3)]
        records = []
        for family,ti,vi in subsets:
            a = [v for i,(_,v) in enumerate(arms["treated"]) if i!=ti]
            b = [v for i,(_,v) in enumerate(arms["vehicle"]) if i!=vi]
            d,p,e,n = exact_test(a,b)
            record = dict(outcome=name,sensitivity=family,
                deleted_treated_cell="" if ti is None else arms["treated"][ti][0],
                deleted_vehicle_cell="" if vi is None else arms["vehicle"][vi][0],
                n_treated=len(a),n_vehicle=len(b),mean_treated=fmean(a),mean_vehicle=fmean(b),
                treated_minus_vehicle=d,percent_lower=100*(1-fmean(a)/fmean(b)),
                exact_p_two_sided=p,extreme_allocations=e,allocations=n)
            rows.append(record)
            records.append(record)
        full = records[0]
        old = next(x for x in baseline["contrasts"] if x["outcome"]==name)
        assert math.isclose(full["treated_minus_vehicle"],old["treated_minus_vehicle"],abs_tol=1e-12)
        assert full["exact_p_two_sided"]==old["exact_permutation_p_two_sided"]
        envelopes = []
        for family in ("delete_one_single_arm","delete_one_each_arm"):
            rs = [r for r in records if r["sensitivity"]==family]
            envelopes.append(dict(sensitivity=family,subset_count=len(rs),negative_difference_count=sum(r["treated_minus_vehicle"]<0 for r in rs),
                difference_range=[min(r["treated_minus_vehicle"] for r in rs),max(r["treated_minus_vehicle"] for r in rs)],
                percent_lower_range=[min(r["percent_lower"] for r in rs),max(r["percent_lower"] for r in rs)],
                exact_p_range=[min(r["exact_p_two_sided"] for r in rs),max(r["exact_p_two_sided"] for r in rs)]))
        neuronal.append(dict(outcome=name,full=full,deletion_envelopes=envelopes,
            all_nine_cross_arm_pairs_negative=all(t[1]<v[1] for t in arms["treated"] for v in arms["vehicle"]),
            minimum_vehicle_minus_maximum_treated=min(v[1] for v in arms["vehicle"])-max(t[1] for t in arms["treated"])))
    book.close()
    running = 0.0
    for rank,i in enumerate(sorted(range(2),key=lambda i:neuronal[i]["full"]["exact_p_two_sided"])):
        running = max(running,min(1.0,(2-rank)*neuronal[i]["full"]["exact_p_two_sided"]))
        neuronal[i]["exact_p_holm_two_outcomes"] = running
    pk = json.loads(paths["exposure_config"].read_text())
    old = json.loads(paths["exposure_summary"].read_text())
    assert old["config_sha256"]==digest(paths["exposure_config"])
    method = design["exposure_methods"]
    nominal = method["nominal_culture_mM"]
    assert nominal==cfg["treatment"]["concentration_uM"]/1000
    parameters = list(itertools.product(pk["molecular_weight_values_mg_per_mmol"],pk["adult_single_dose_cmax_ug_ml"],pk["food_peak_factors"],pk["human_plasma_unbound_fraction_bounds"],pk["adult_half_life_hours"]))
    assert len(parameters)==64
    exposure_rows, thresholds = [], []
    for time in method["post_peak_h"]:
        selected = []
        for mw,peak,food,fu,half in parameters:
            free = peak*pk["food_peak_factors"][food]/mw*fu*2**(-time/half)
            row = dict(post_peak_h=time,mass_convention=mw,reported_adult_peak_ug_ml=peak,food_context=food,plasma_fu=fu,half_life_h=half,modeled_free_mM=free,maximum_culture_fu_for_overlap=min(1.0,free/nominal))
            selected.append(row)
            exposure_rows.append(row)
        thresholds.append(dict(post_peak_h=time,modeled_free_mM_range=[min(r["modeled_free_mM"] for r in selected),max(r["modeled_free_mM"] for r in selected)],
            maximum_culture_fu_requirement_range=[min(r["maximum_culture_fu_for_overlap"] for r in selected),max(r["maximum_culture_fu_for_overlap"] for r in selected)],
            selected_conditions=[dict(culture_fu=x,overlap_in_all_assumptions=all(r["modeled_free_mM"]>=nominal*x for r in selected),overlap_in_any_assumption=any(r["modeled_free_mM"]>=nominal*x for r in selected)) for x in method["culture_fu"]]))
    assert all(math.isclose(a,b,abs_tol=1e-12) for a,b in zip(thresholds[0]["modeled_free_mM_range"],old["adult_free_peak_mM_bounded_range"]))
    tails = []
    for x in method["culture_fu"]:
        times = [max(0.0,half*math.log2(peak*pk["food_peak_factors"][food]/mw*fu/(nominal*x))) for mw,peak,food,fu,half in parameters]
        tails.append(dict(culture_fu=x,post_peak_overlap_h_range=[min(times),max(times)]))
    out = root/args.out
    out.mkdir(parents=True,exist_ok=True)
    write_tsv(out/"deletion_results.tsv",rows)
    write_tsv(out/"exposure_thresholds.tsv",exposure_rows)
    result = dict(schema_version=1,status="exploratory_public_sensitivity_complete",independent_new_biological_units=0,
        engine_sha256=digest(Path(__file__)),design_sha256=digest(design_path),inputs=design["inputs"],neuronal=neuronal,
        exposure=dict(nominal_culture_mM=nominal,assumptions_per_time=64,post_peak_requirements=thresholds,duration_envelopes=tails,model_scope=method["scope"]),
        unresolved_endpoints=design["unresolved_endpoints"])
    (out/"summary.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(dict(out=str(out),summary_sha256=digest(out/"summary.json"),deletion_rows=len(rows),threshold_rows=len(exposure_rows))))


if __name__=="__main__":
    main()
