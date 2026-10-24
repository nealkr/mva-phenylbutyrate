#!/usr/bin/env python3
"""Reanalyse fixed public mouse construct/stressor chromosome-alignment data."""
import argparse
import csv
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
import openpyxl
from scipy.stats import t


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def bh(values):
    p = np.asarray(values, float)
    order = np.argsort(p)
    adjusted = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    out = np.empty(len(p))
    out[order] = np.minimum(adjusted, 1)
    return out.tolist()


def exact_null(x, y):
    pooled = np.r_[x, y]
    observed = float(x.mean() - y.mean())
    values = []
    for allocation in itertools.combinations(range(len(pooled)), len(x)):
        select = np.zeros(len(pooled), dtype=bool)
        select[list(allocation)] = True
        values.append(float(pooled[select].mean() - pooled[~select].mean()))
    return {"allocations": len(values), "two_sided_p": sum(abs(v) >= abs(observed) - 1e-12
             for v in values) / len(values), "null": "exchangeable_line_outcomes_across_construct_groups"}


def bootstrap_diff(x, y, rng, draws):
    xm = rng.choice(x, size=(draws, len(x))).mean(axis=1)
    ym = rng.choice(y, size=(draws, len(y))).mean(axis=1)
    return np.quantile(xm - ym, [0.025, 0.975]).tolist()


def pairing_envelope(vehicle, stressor):
    """Conditional on reuse of six lines, enumerate every missing line pairing."""
    critical = float(t.ppf(0.975, len(vehicle) - 1))
    observed = float(stressor.mean() - vehicle.mean())
    errors = []
    for order in itertools.permutations(range(len(vehicle))):
        delta = stressor - vehicle[list(order)]
        errors.append(float(delta.std(ddof=1) / np.sqrt(len(delta))))
    se_min, se_max = min(errors), max(errors)
    return {"assumption": "same_six_lines_in_each_arm_with_unknown_correspondence",
            "pairings_enumerated": len(errors), "mean_change_pp": observed,
            "paired_t_ci95_outer_envelope_pp": [observed - critical * se_max, observed + critical * se_max],
            "paired_standard_error_range_pp": [se_min, se_max]}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--design", default="data/public/bubr1-functional-v2/fidelity-config.json")
    parser.add_argument("--output", default="results/bubr1_fidelity_calibration")
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args()
    project = args.project_root.resolve()
    config_path = project / args.design
    config = json.loads(config_path.read_text())
    source = project / config["source"]["path"]
    assert sha(source) == config["source"]["sha256"], "Frozen source mismatch"
    workbook = openpyxl.load_workbook(source, read_only=True, data_only=True)
    sheet = workbook[config["source"]["sheet"]]
    rows = list(sheet.iter_rows(values_only=True))
    names = config["constructs"]
    assert list(rows[3][2:7]) == names and list(rows[15][2:7]) == names
    assert "Vehicle" in rows[3][1] and "10 nM AZD" in rows[15][1]
    vehicle = np.asarray([row[2:7] for row in rows[4:10]], float).T
    stressor = np.asarray([row[2:7] for row in rows[16:22]], float).T
    assert vehicle.shape == stressor.shape == (5, 6)
    assert np.isfinite(vehicle).all() and np.isfinite(stressor).all()
    assert ((vehicle >= 0) & (vehicle <= 100)).all() and ((stressor >= 0) & (stressor <= 100)).all()
    assert np.allclose(vehicle.mean(axis=1), np.asarray(rows[10][2:7], float))
    assert np.allclose(stressor.mean(axis=1), np.asarray(rows[22][2:7], float))
    assert config["biological_unit"]["n_per_construct_arm"] == 6
    assert config["biological_unit"]["status"] == "six_independently_generated_mef_lines_supported_by_figure_legend_and_methods"
    if args.verify_only:
        print(json.dumps({"verified": True, "source_sha256": sha(source), "line_observations": 60,
                          "stressor_means_pp": dict(zip(names, stressor.mean(axis=1).tolist()))}))
        return
    rng = np.random.default_rng(config["seed"])
    primary = []
    for index in range(1, len(names)):
        null = exact_null(stressor[index], stressor[0])
        primary.append({"construct": names[index], "reference": names[0],
                        "effect_pp": float(stressor[index].mean() - stressor[0].mean()),
                        "ci95_pp": bootstrap_diff(stressor[index], stressor[0], rng, config["bootstrap_draws"]),
                        **null})
    for row, q in zip(primary, bh([row["two_sided_p"] for row in primary])):
        row["bh_q_four_primary_comparisons"] = q
    groups = []
    for name, v, s in zip(names, vehicle, stressor):
        groups.append({"construct": name, "vehicle_line_values_pp": v.tolist(),
                       "stressor_line_values_pp": s.tolist(), "vehicle_mean_pp": float(v.mean()),
                       "stressor_mean_pp": float(s.mean()),
                       "independent_arm_change_ci95_pp": bootstrap_diff(s, v, rng, config["bootstrap_draws"]),
                       "unknown_pairing_sensitivity": pairing_envelope(v, s)})
    interactions = []
    wt = groups[0]["unknown_pairing_sensitivity"]
    for row in groups[1:]:
        env = row["unknown_pairing_sensitivity"]
        effect = env["mean_change_pp"] - wt["mean_change_pp"]
        # df5 critical is conservative for two independent n6 genotype line groups,
        # over every possible within-group vehicle/stressor correspondence.
        margin = float(t.ppf(0.975, 5)) * np.hypot(env["paired_standard_error_range_pp"][1],
                                                 wt["paired_standard_error_range_pp"][1])
        interactions.append({"construct": row["construct"], "difference_in_drug_changes_vs_wt_pp": effect,
                             "conditional_unknown_pairing_ci95_outer_envelope_pp": [effect-margin, effect+margin],
                             "status": "descriptive_sensitivity_not_primary_test"})
    result = {"schema_version": 1, "design_sha256": sha(config_path), "engine_sha256": sha(Path(__file__)),
              "source_sha256": sha(source), "species": "Mus musculus",
              "outcome": "MEF-line percentage of cells with one or more misaligned chromosomes after monastrol washout",
              "primary": primary, "groups": groups, "interactions": interactions,
              "bootstrap_draws": config["bootstrap_draws"],
              "canonical_design_sha256": config.get("canonical_design_sha256"),
              "interpretation": config.get("interpretation", []), "uncertainty": config.get("uncertainty", [])}
    out = project / args.output
    out.mkdir(parents=True, exist_ok=True)
    (out / "summary.json").write_text(json.dumps(result, indent=2) + "\n")
    with (out / "line_outcomes.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["construct", "arm", "source_row", "misalignment_percent", "biological_unit"])
        for index, name in enumerate(names):
            for arm, array, start in [("vehicle", vehicle, 5), ("AZD1152_10nM", stressor, 17)]:
                for i, value in enumerate(array[index]):
                    writer.writerow([name, arm, start+i, value, "MEF_line_mean_no_deposited_line_ID"])
    print(json.dumps({"output": str(out), "primary": primary}, separators=(",", ":")))


if __name__ == "__main__":
    main()
