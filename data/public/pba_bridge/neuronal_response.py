#!/usr/bin/env python3
"""Reanalyze two public neuronal response panels using differentiation units."""
import argparse
import hashlib
import json
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--config', type=Path, default=Path(__file__).with_name('neuronal_response_config.json'))
    ap.add_argument('--source', type=Path)
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    cfg = json.loads(args.config.read_text())
    source = args.source or args.config.with_name(cfg['source_file'])
    source_hash = hashlib.sha256(source.read_bytes()).hexdigest()
    if source_hash != cfg['source_sha256']:
        raise ValueError('Public source workbook hash mismatch')
    data = pd.read_excel(source, sheet_name=cfg['sheet'], header=None)
    observations, contrasts = [], []
    for outcome in cfg['outcomes']:
        cell = outcome['panel_cell']
        if str(data.iloc[int(cell[1:]) - 1, ord(cell[0]) - ord('A')]) != outcome['panel']:
            raise ValueError(f"Panel row-map mismatch: {outcome['id']}")
        vals = {}
        for group, meta in cfg['groups'].items():
            column = ord(meta['column']) - ord('A')
            values = []
            for i, row in enumerate(outcome['data_rows'], 1):
                value = float(data.iloc[row - 1, column])
                if not np.isfinite(value) or value < 0:
                    raise ValueError('Nonfinite or negative source intensity')
                values.append(value)
                observations.append({'outcome': outcome['id'], 'group': group,
                                     'source_cell': f"{meta['column']}{row}",
                                     'source_position': i, 'value': value})
            if len(values) != cfg['replicates_per_arm']:
                raise ValueError('Unexpected differentiation count')
            vals[group] = np.array(values)
        a, b = vals['treated'], vals['vehicle']
        va, vb = a.var(ddof=1) / len(a), b.var(ddof=1) / len(b)
        se = np.sqrt(va + vb)
        df = (va + vb) ** 2 / (va ** 2 / (len(a) - 1) + vb ** 2 / (len(b) - 1))
        delta = a.mean() - b.mean()
        radius = stats.t.ppf(0.975, df) * se
        test = stats.ttest_ind(a, b, equal_var=False)
        pooled = np.concatenate([a, b])
        permutation_deltas = []
        for allocation in combinations(range(len(pooled)), len(a)):
            mask = np.zeros(len(pooled), dtype=bool)
            mask[list(allocation)] = True
            permutation_deltas.append(pooled[mask].mean() - pooled[~mask].mean())
        extreme = int(sum(abs(x) >= abs(delta) - 1e-12 for x in permutation_deltas))
        exact_p = extreme / len(permutation_deltas)
        contrasts.append({'outcome': outcome['id'], 'n_treated': len(a), 'n_vehicle': len(b),
                          'mean_reference': vals['reference'].mean(), 'mean_vehicle': b.mean(),
                          'mean_treated': a.mean(), 'treated_minus_vehicle': delta,
                          'difference_ci_low': delta - radius, 'difference_ci_high': delta + radius,
                          'percent_lower_than_vehicle': 100 * (1 - a.mean() / b.mean()),
                          'treated_minus_reference_descriptive': a.mean() - vals['reference'].mean(),
                          'welch_df': df, 'p_nominal': float(test.pvalue),
                          'exact_permutation_p_two_sided': exact_p,
                          'exact_permutation_allocations': len(permutation_deltas),
                          'exact_permutation_extreme_allocations': extreme,
                          'minimum_two_sided_permutation_p': 2 / len(permutation_deltas)})
    for input_key, output_key in [('p_nominal', 'p_holm_two_outcomes'),
                                  ('exact_permutation_p_two_sided', 'exact_permutation_p_holm_two_outcomes')]:
        order = np.argsort([r[input_key] for r in contrasts])
        running = 0.0
        for rank, index in enumerate(order):
            running = max(running, min(1.0, (len(order) - rank) * contrasts[index][input_key]))
            contrasts[index][output_key] = running
    args.out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(observations).to_csv(args.out / 'observations.tsv', sep='\t', index=False)
    pd.DataFrame(contrasts).to_csv(args.out / 'contrasts.tsv', sep='\t', index=False)
    summary = {'analysis': 'public source-paper response reanalysis, exploratory',
               'source_sha256': source_hash,
               'engine_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               'config_sha256': hashlib.sha256(args.config.read_bytes()).hexdigest(),
               'observation_count': len(observations), 'comparison_count': len(contrasts),
               'biological_unit': cfg['unit'], 'analysis_plan': cfg['analysis'],
               'groups': cfg['groups'], 'contrasts': contrasts}
    (args.out / 'summary.json').write_text(json.dumps(summary, indent=2))
    print(json.dumps({'out': str(args.out), 'n_observations': len(observations),
                      'results': [{'outcome': r['outcome'], 'percent_lower': round(r['percent_lower_than_vehicle'], 3),
                                   'difference_ci': [round(r['difference_ci_low'], 4), round(r['difference_ci_high'], 4)],
                                   'p_holm': round(r['p_holm_two_outcomes'], 6),
                                   'exact_permutation_p': r['exact_permutation_p_two_sided'],
                                   'exact_permutation_allocations': r['exact_permutation_allocations']} for r in contrasts]}))


if __name__ == '__main__':
    main()
