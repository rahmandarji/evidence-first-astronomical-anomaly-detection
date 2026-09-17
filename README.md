# Evidence-First Astronomical Anomaly Detection

Evidence-first, multi-view research framework for prioritizing unusual astronomical sources while separating statistical novelty, evidence quality, systematic risk, and astrophysical interpretation.

> Unusual does not mean novel, real, or confirmed.

## Research project

This repository contains the public research code and selected notebooks for an independent astronomy research project using NEOWISE single-exposure data.

The project asks:

> Can a multi-view, uncertainty-aware novelty detection system prioritize astronomical sources whose observed behavior is robustly unusual and poorly explained by known source classes or tested observational/systematic effects?

The research follows an evidence-first pipeline:

```text
UNUSUAL
  ↓
DATA-QUALITY VERIFIED
  ↓
ARTIFACT / BLENDING / MOTION CHECKS
  ↓
CROSS-SURVEY / CATALOG CHECKS
  ↓
KNOWN-CLASS CHECKS
  ↓
ASTROPHYSICALLY PLAUSIBLE
  ↓
FOLLOW-UP CANDIDATE
  ↓
OBSERVATIONAL CONFIRMATION
A candidate being unusual is not itself evidence of a new astronomical class or new physics.
Data source
The primary data source is the NEOWISE Single-exposure Source Database distributed through the NASA/IPAC Infrared Science Archive (IRSA).
NEOWISE provides repeated W1 (3.4 μm) and W2 (4.6 μm) measurements together with positions, uncertainties, observation times, and quality information.
Official data documentation is available through NASA/IPAC IRSA.
The public repository does not contain the full NEOWISE survey archive. The survey data are extremely large, so users should obtain only the data required for their own reproduction.
Research field
The development analysis used a small sky region centered approximately on:
RA  ≈ 179.9–180.1 degrees
Dec ≈ -0.1–+0.1 degrees
The development field should not be interpreted as a survey-wide statistical sample.
Expected local data layout
Large survey datasets are intentionally kept outside Git.
For local reproduction, data may be organized as:
data/
├── raw/
├── processed/
└── results/
A project-specific input may look like:
data/
└── raw/
    └── exp005_region_detections.csv
The exact data products required by individual historical experiments may differ.
Do not commit raw NEOWISE survey files, large Parquet datasets, downloaded image archives, or the private research archive.
Environment
The required Python packages are listed in:
requirements.txt
Install them with:
python -m pip install -r requirements.txt
Core dependencies include:
NumPy
Pandas
SciPy
Matplotlib
Astropy
PyArrow
scikit-learn
Requests
Public notebook order
The public notebooks are a curated subset of the larger private research history.
01 — Initial NEOWISE validation
notebooks/public/01_initial_neowise_validation.ipynb
Early validation of the NEOWISE candidate-selection and photometric workflow.
02 — Validation population
notebooks/public/02_validation_population.ipynb
Construction and validation of the source population used for later comparisons.
03 — Image analysis
notebooks/public/03_image_analysis.ipynb
Image-level investigation of candidate observations, source morphology, photometry, and artifact-related evidence.
04 — Phase IV ML analysis
notebooks/public/04_phase4_ml_analysis.ipynb
Feature construction and exploratory unsupervised novelty analysis.
The machine-learning analysis is used for triage and diagnosis. It is not presented as confirmation of astrophysical novelty.
05 — Analysis freeze audit
notebooks/public/05_analysis_freeze_audit.ipynb
Audit and provenance checks associated with the frozen analysis state.
Reusable source code
Reusable utilities are provided under:
src/
├── coordinates.py
├── statistics.py
└── validation.py
These modules provide coordinate utilities, uncertainty-aware statistical calculations, and variability-validation helpers extracted from the broader research workflow.
Analysis architecture
The project separates three major concepts:
Novelty
Evidence quality
Systematic risk
These are not collapsed into a single arbitrary weighted score.
The broader framework uses multiple evidence views:
temporal behavior
photometric behavior
image evidence
observational/systematic evidence
external survey information
known-class similarity
Missing information is represented as unknown rather than automatically interpreted as normal or negative evidence.
Analysis status
The research analysis was frozen through EXP-056.
The final candidate set was:
Source
Priority
Interpretation
1978
Tier 1 — highest
Persistent unusual behavior; highest follow-up priority
1974
Tier 1 — high
Persistent unusual behavior; high follow-up priority
233
Tier 2
Requires additional controls
1976
Tier 2
Requires additional controls
234
Control/sanity
Control-compatible
1979
Control/sanity
Control-compatible
No candidate is confirmed as a new astronomical class or new physical phenomenon.
The final analysis identifies 1978 and 1974 as the highest-priority sources for additional investigation and follow-up.
Important candidate evidence
Source 1978
Source 1978 showed persistent excess W1 variability relative to the uncertainty-aware baseline.
The Phase III validator gave a W1 reduced chi-square of approximately 3.94, with the later Phase VII recomputation giving approximately 4.14 under a slightly different direct reconstruction.
Leave-one-out testing retained the signal, with a minimum reduced chi-square of approximately 2.80.
Matched-control comparisons found the source unusual within the tested control population.
A nearby NEOWISE source was identified, so blending or contamination remains a live hypothesis. Direct neighbor tracking did not provide a strong reproduction of the target behavior.
No valid SIMBAD identification was returned, while the Gaia query was unavailable. These are unknown states, not evidence that the source has no identity.
The source is therefore a high-priority unexplained candidate, not a confirmed astrophysical discovery.
Source 1974
Source 1974 also showed persistent excess W1 variability.
The Phase III validator gave a W1 reduced chi-square of approximately 3.03, with the later Phase VII recomputation giving approximately 3.05.
Leave-one-out testing retained the signal, with a minimum reduced chi-square of approximately 2.14.
Matched-control comparisons found the source unusual within the tested control population.
A nearby NEOWISE source was identified, so blending remains a possible explanation, although direct neighbor tracking did not strongly reproduce the target behavior.
No valid SIMBAD identification was returned and the Gaia query was unavailable.
The source is therefore a high-priority unexplained candidate, not a confirmed astrophysical discovery.
Scientific interpretation
The analysis supports prioritizing unusual sources for further investigation.
It does not establish:
a new astronomical class;
new physics;
confirmed intrinsic variability;
absence of blending;
absence of instrumental/systematic effects;
a formal survey-wide false-positive rate;
observational confirmation.
The statement "unexplained" means that the tested explanations did not adequately account for the observations available in this analysis. It does not mean that all possible explanations have been eliminated.
Public/private boundary
This repository is the public research-code and reproducibility layer.
The complete private research archive remains outside this repository.
The private archive contains the broader experimental history, intermediate artifacts, large data products, frozen evidence archive, and working materials that are not required for the public code release.
The public repository therefore does not claim that a fresh clone reconstructs every historical experiment byte-for-byte.
The public release provides:
reusable research code;
selected research notebooks;
dependency information;
data-source documentation;
methodology;
provenance information;
selected analysis context.
What is not included
The following are intentionally excluded from the public repository:
full NEOWISE raw/survey datasets;
large downloaded Parquet archives;
downloaded image archives;
private Google Drive paths;
private research archives;
personal working files;
credentials or secrets;
uncurated intermediate experiment outputs.
Reproducibility expectations
A researcher cloning this repository should treat the notebooks as research records and methodology references as well as executable code.
Some historical notebooks depend on data products generated during earlier experiments. Those data products are not all included in this public repository.
Reproduction of a specific historical result may therefore require obtaining the corresponding public survey data and reconstructing the relevant intermediate data products.
The frozen numerical conclusions documented in the research archive should not be silently regenerated or replaced with new values merely to make the public repository appear more complete.
Scientific limitations
Important limitations include:
the development field is small and does not constitute a survey-wide population study;
the validator-positive set is not a formal false-positive rate;
matched-control comparisons are limited in size and field;
the control population is not equivalent to a statistically complete survey population;
missing external information is not treated as negative evidence;
image-level source-like structure does not prove astrophysical variability;
neighboring-source tests do not completely eliminate blending;
unsupervised ML novelty is not equivalent to astrophysical novelty;
exact ML ranks can be sensitive to population and feature sampling;
no candidate has observational confirmation in this analysis.
Evidence language
The project uses conservative evidence labels.
Preferred language includes:
CONFIRMED
SUPPORTED
LIKELY
PLAUSIBLE
INCONCLUSIVE
REJECTED
NOT TESTED
UNKNOWN
These labels describe the strength or availability of evidence and should not be interpreted as probabilities unless a statistical probability model has explicitly been constructed.
Citation
If you use the research code or methodology, please cite this repository.
Citation metadata are provided in:
CITATION.cff
License
This project is released under the MIT License.
See CITATION.cff for author and citation information.
Research principle
Detect → Verify → Explain
An anomaly is a hypothesis, not a discovery.
