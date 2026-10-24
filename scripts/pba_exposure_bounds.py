#!/usr/bin/env python3
"""Compute explicit free-exposure requirements under bounded public PK inputs."""
import argparse
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    config_path = Path(args.config)
    config = json.loads(config_path.read_text())
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    for key in ('human_plasma_unbound_fraction_bounds', 'culture_unbound_fraction_sensitivity'):
        if not all(0 < x <= 1 for x in config[key]):
            raise ValueError(f'Invalid fractions: {key}')
    rows = []
    for mw, peak, half, food, plasma_fu, culture_fu, nominal in itertools.product(
            config['molecular_weight_values_mg_per_mmol'],
            config['adult_single_dose_cmax_ug_ml'], config['adult_half_life_hours'],
            config['food_peak_factors'], config['human_plasma_unbound_fraction_bounds'],
            config['culture_unbound_fraction_sensitivity'], config['culture_nominal_mM']):
        total_peak = peak * config['food_peak_factors'][food] / mw
        free_peak = total_peak * plasma_fu
        assumed_free_culture = nominal * culture_fu
        tail_hours = (half * math.log2(free_peak / assumed_free_culture)
                      if free_peak > assumed_free_culture else 0.0)
        rows.append({'molecular_weight_mg_per_mmol_convention': mw,
                     'adult_fasted_reported_mean_peak_ug_ml': peak, 'food_context': food,
                     'total_peak_mM': total_peak, 'plasma_fu_bound': plasma_fu,
                     'free_peak_mM_bound': free_peak, 'half_life_h_bound': half,
                     'culture_nominal_mM': nominal, 'culture_fu_assumption': culture_fu,
                     'culture_free_mM_assumption': assumed_free_culture,
                     'maximum_culture_fu_for_peak_overlap': min(1.0, free_peak / nominal),
                     'post_peak_tail_above_assumed_culture_free_h': tail_hours,
                     'peak_overlap_under_assumptions': free_peak >= assumed_free_culture})
    with (out / 'exposure_scenarios.tsv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    summary = {'config_sha256': hashlib.sha256(config_path.read_bytes()).hexdigest(),
               'scenario_rows': len(rows), 'independent_experiments': 0,
               'adult_total_peak_mM_range': [min(x['total_peak_mM'] for x in rows), max(x['total_peak_mM'] for x in rows)],
               'adult_free_peak_mM_bounded_range': [min(x['free_peak_mM_bound'] for x in rows), max(x['free_peak_mM_bound'] for x in rows)],
               'culture_binding_dependence': [], 'scope': config['model'],
               'interpretation': config['inference'],
               'pk_range_scope': config['pk_range_scope'],
               'mass_conversion_scope': config['mass_conversion_basis']}
    for nominal in config['culture_nominal_mM']:
        selected = [x for x in rows if x['culture_nominal_mM'] == nominal]
        summary['culture_binding_dependence'].append({
            'nominal_mM': nominal,
            'maximum_allowable_culture_fu_for_peak_overlap_range': [min(x['maximum_culture_fu_for_peak_overlap'] for x in selected), max(x['maximum_culture_fu_for_peak_overlap'] for x in selected)],
            'assumptions_with_peak_overlap': sum(x['peak_overlap_under_assumptions'] for x in selected),
            'sensitivity_rows': len(selected),
            'maximum_post_peak_tail_h_across_assumptions': max(x['post_peak_tail_above_assumed_culture_free_h'] for x in selected)})
    (out / 'summary.json').write_text(json.dumps(summary, indent=2))
    print(json.dumps({'scenario_rows': len(rows), 'total_peak_range': summary['adult_total_peak_mM_range'], 'free_peak_bounds': summary['adult_free_peak_mM_bounded_range']}))


if __name__ == '__main__':
    main()
