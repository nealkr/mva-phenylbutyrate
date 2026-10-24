#!/usr/bin/env python3
"""Diagnostic ranking comparison on fixed, public BUBR1 experimental controls."""
import argparse
import csv
import hashlib
import json
from pathlib import Path


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def auc(positive, negative):
    comparisons = [(x > y) + 0.5*(x == y) for x in positive for y in negative]
    return {"auroc": sum(comparisons)/len(comparisons), "concordant_pair_equivalents": sum(comparisons),
            "positive_negative_pairs": len(comparisons), "positive_count": len(positive),
            "negative_count": len(negative), "higher_score_direction": "reduced_abundance"}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--config", default="data/public/bubr1-functional-v2/prediction-config.json")
    parser.add_argument("--output", default="results/bubr1_prediction_calibration")
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args()
    project = args.project_root.resolve()
    config_path = project / args.config
    config = json.loads(config_path.read_text())
    input_path = project / config["input"]
    assert sha(input_path) == config["input_sha256"]
    with input_path.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["variant"] for row in rows] == config["control_variants"]
    assert len(rows) == len({row["variant"] for row in rows}) == 14
    result = []
    for index, expanded in enumerate([False, True]):
        selected = rows if expanded else [row for row in rows if int(row["expanded_only"]) == 0]
        positive = [row for row in selected if int(row["reduced_abundance_label"]) == 1]
        negative = [row for row in selected if int(row["reduced_abundance_label"]) == 0]
        assert len(selected) == config["expected_control_counts"][index]
        assert [len(positive), len(negative)] == config["expected_class_counts"][index]
        scores = {column: auc([float(row[column]) for row in positive],
                              [float(row[column]) for row in negative]) for column in config["predictors"]}
        result.append({"control_set": "expanded14" if expanded else "initial12",
                       "predictors": scores, "unique_residue_count": len({row["residue"] for row in selected}),
                       "positive_unique_residues": len({row["residue"] for row in positive}),
                       "negative_unique_residues": len({row["residue"] for row in negative}),
                       "variants": [row["variant"] for row in selected]})
    output = {"schema_version": 1, "engine_sha256": sha(Path(__file__)), "config_sha256": sha(config_path),
              "input_sha256": sha(input_path), "parent_prediction_sha256": config["parent_prediction_sha256"],
              "comparisons": result,
              "data_type": "previously_published_experimental_controls_and_precomputed_prediction_scores",
              "new_model_predictions": False, "bootstrap_or_p_values": False,
              "unit": "variant_controls_with_shared_residue_clusters",
              "interpretation": "diagnostic_rankings_for_heterogeneous_cellular_abundance_labels"}
    if not args.verify_only:
        out = project / args.output
        out.mkdir(parents=True, exist_ok=True)
        (out / "summary.json").write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, separators=(",", ":")))


if __name__ == "__main__":
    main()
