# MVA drug-repurposing analyses

Neal Krishna | Signal Works

Eight public-data analyses support the research proposal. Two optional analyses examine mitochondrial import and sensitivity to individual observations or exposure assumptions.

## Run

Python 3.10 was used with the dependencies pinned in `requirements.txt`.

```sh
python -m venv .venv
```

Activate the environment with `. .venv/bin/activate` on Linux/macOS or `.venv\Scripts\Activate.ps1` in Windows PowerShell, then run:

```sh
python -m pip install -r requirements.txt
python fetch_inputs.py --group core
python run.py --core --check
python verify.py
```

Run commands from the repository root. Source downloads require internet access.

| Analysis | Comparison |
| --- | --- |
| `pba_exposure` | Phenylbutyrate exposure under explicit mass, binding, and food-effect assumptions. |
| `bubr1_fidelity` | BUBR1 construct-by-stressor chromosome alignment. |
| `prediction_controls` | AlphaMissense and RaSP rankings against published BUBR1 abundance controls. |
| `pba_neuronal` | Two neuronal endpoints comparing 100 micromolar phenylbutyrate with vehicle. |
| `niacin_cachexia` | Mouse muscle function and energy endpoints. |
| `nac_fly` | Respiration and ROS across three fly haplotypes. |
| `nac_ndufa3` | Chromosome-error endpoints in the NDUFA3 experiment. |
| `nac_g6pd` | Chromosome-fidelity comparisons in G6PD/BAG3 experiments. |

The 12- and 14-control sets overlap. The neuronal endpoints share one source experiment; each has three differentiation summaries per arm.

## Optional analyses

```sh
python fetch_inputs.py --group import
python run.py --analysis import_kinetics --check
python run.py --analysis submission_sensitivity --check
```

The import analysis reproduces source-paired mitochondrial trajectories. The sensitivity analysis reuses the core inputs for 32 deletion records and 320 exposure scenarios.

## Reproducibility

- `provenance.json` records source URLs, input hashes, rights, and preparation steps.
- `analyses.json` lists the scripts, configurations, and inputs for each analysis.
- `expected/` holds reference results. The `--check` option compares new outputs in `results/` with these references.
- `manifest.json` records release-file hashes. `verify.py` checks integrity.

## License and citation

Original materials use [CC BY 4.0](LICENSE). Third-party data retain their own licences; some inputs are fetched separately. Cite this package using [CITATION.cff](CITATION.cff) and the original sources in [THIRD_PARTY.md](THIRD_PARTY.md).

## Acknowledgment

This work was made possible through the Hackathon, organized by Sage Bionetworks in partnership with the MVA Society, Hugging Face, and BEACON (The Benchmarking, Evaluation, and Assessment Consortium for Science), with prize sponsorship from AWS and Anthropic. We are deeply grateful to the child and their family who generously contributed their data and their story to advance research into this rare disease. We acknowledge their trust in making this Hackathon possible.
