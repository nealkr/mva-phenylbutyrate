# Third-party sources and rights

The package's CC BY 4.0 license covers original material by Neal Krishna / Signal Works. It does not relicense third-party sources. `provenance.json` records the exact source URLs, hashes, included paths, and transformations; the notices below explain how those sources are handled.

## Scientific and regulatory data

| Source | Rights and package treatment |
| --- | --- |
| [Weaver et al., eLife 5:e16620 (2016), Figure 7](https://doi.org/10.7554/eLife.16620) | CC BY 4.0. Figure 7 source data support the construct/stressor comparison; extraction and analysis are identified in provenance. |
| [Beltrà et al., Nature Communications (2023)](https://doi.org/10.1038/s41467-023-37595-6) | CC BY 4.0. The source workbook is fetched separately; selected numerical observations support the niacin comparisons. |
| [Rodriguez et al., UCL source dataset, version 1](https://doi.org/10.5522/04/16539723.v1) | CC0. Source data support the fly NAC respiration/ROS comparisons across WT/COX/BAR haplotypes; each independent chamber contains two larvae. Provenance identifies the selected measurements. |
| [Wang et al., EMBO Reports (2026)](https://doi.org/10.1038/s44319-026-00937-8) | Source data are CC0. Selected data support the NDUFA3/NAC comparison. Article text and figures are not redistributed. |
| [Cheng et al., Nature Metabolism (2023)](https://doi.org/10.1038/s42255-023-00832-9) | Source workbooks are under exclusive publisher rights and are fetched from the publisher rather than bundled. The fetch/preparation step extracts selected numerical observations to CSV, with source locations and derivation recorded in provenance. These generated CSVs are not bundled. The CC BY 4.0 notice does not apply to the fetched workbooks. |
| [Nambara et al., Communications Biology (2026)](https://doi.org/10.1038/s42003-026-09902-6) | Article license: CC BY-NC-ND 4.0. The source workbook is fetched rather than bundled. The analysis records direct numerical observations with source locations for the two 100-micromolar PBA panels. The article, workbook, and figures are not redistributed; the original-material license does not remove source restrictions. |
| [Amponsah et al., Nature Communications (2025)](https://doi.org/10.1038/s41467-025-60857-4) | CC BY 4.0. The source workbook for the optional import-kinetics replication is fetched separately. |
| [PHEBURANE prescribing information](https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=d15f2cbf-4d22-40b5-8ef9-16b9d323a7f9) and [FDA parent-PBA binding information](https://www.accessdata.fda.gov/drugsatfda_docs/label/2013/203284s000lbl.pdf) | Public regulatory sources for numerical exposure assumptions. The package records numerical observations and source links; the original label documents are not bundled or relicensed. |

Numerical extraction and statistical reanalysis are performed for the comparisons identified in `analyses.json`. The resulting calculations are separate from source figures and workbooks. Related panels from one paper are not represented as independent experimental replication.

## Prediction controls

The bundled label table records selected factual abundance observations from [Suijkerbuijk et al. (2010)](https://doi.org/10.1158/0008-5472.CAN-09-4319), [Suijkerbuijk et al. (2012)](https://doi.org/10.1016/j.devcel.2012.03.009), and [Gama Braga et al. (2020)](https://doi.org/10.1016/j.celrep.2020.108397). Article text and figures are not bundled. These heterogeneous labels are literature-selected abundance controls, not a prospective benchmark of protein stability.

[AlphaMissense precomputed predictions](https://github.com/google-deepmind/alphamissense) used in the control table are CC BY 4.0. Copyright 2023 DeepMind Technologies Limited. The bundled material is an extracted subset with unchanged score values; source attribution and the extraction record are retained in `provenance.json`.

[RaSP](https://doi.org/10.7554/eLife.82593) precomputed scores come from the authors' [ERDA data archive](https://sid.erda.dk/sharelink/fFPJWflLeE). A redistribution license for the ERDA score file has not been established. The full CSV, RaSP score column, and source-model pLDDT column are therefore fetch-only. The fetch/preparation step combines those columns with the bundled labels and AlphaMissense subset to reconstruct the retained controls table. The RaSP paper or code license is not treated as a license for the separate data archive.

## Python dependencies

Dependencies are installed separately through `requirements.txt`; their code, wheels, and environments are not bundled.

| Dependency | Version | License and primary notice |
| --- | --- | --- |
| NumPy | 2.2.6 | [BSD 3-Clause](https://raw.githubusercontent.com/numpy/numpy/v2.2.6/LICENSE.txt) |
| SciPy | 1.15.3 | [BSD 3-Clause](https://raw.githubusercontent.com/scipy/scipy/v1.15.3/LICENSE.txt) |
| pandas | 2.3.3 | [BSD 3-Clause](https://raw.githubusercontent.com/pandas-dev/pandas/v2.3.3/LICENSE) |
| openpyxl | 3.1.5 | [MIT](https://pypi.org/project/openpyxl/3.1.5/) |

Any later redistribution of dependency code or binary wheels must retain the applicable upstream notices, including notices for bundled subcomponents.
