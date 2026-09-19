Evidence-First Astronomical Anomaly Detection

«Evidence-first, multi-view research framework for prioritizing unusual astronomical sources while separating statistical novelty, evidence quality, systematic risk, and astrophysical interpretation.»

Research principle:

Detect → Verify → Explain

«An anomaly is a hypothesis, not a discovery.»

---

Overview

This repository contains the public research code, selected notebooks, methodology documentation, and accompanying manuscript for an independent astronomy research project using NEOWISE single-exposure infrared data.

The project investigates whether a multi-view, uncertainty-aware workflow can identify astronomical sources whose observed behavior is unusually difficult to explain by the tested observational, instrumental, catalog, and known-source-class explanations.

The central methodological distinction is:

Unusual ≠ Novel ≠ Real ≠ Confirmed

A statistically unusual measurement is therefore treated as the beginning of an investigation, rather than as evidence of a discovery.

---

📄 Accompanying Paper

The canonical Paper 1 manuscript is included in this repository:

"paper/evidence-First-Astronomical-Anomaly-Investigation.pdf"

Title: Evidence-First Astronomical Anomaly Investigation

The manuscript documents:

- the frozen EXP-056 analysis state;
- the statistical screening methodology;
- robustness and leave-one-out testing;
- image-level investigation;
- observational and systematic checks;
- blending and neighboring-source analysis;
- cross-band behavior;
- external catalog checks;
- candidate evidence dossiers;
- limitations and reproducibility boundaries.

The PDF is the canonical manuscript associated with this research release.

---

Research Question

«Can a multi-view, uncertainty-aware novelty detection system prioritize astronomical sources whose observed behavior is robustly unusual and poorly explained by known source classes or tested observational/systematic effects?»

This question is deliberately narrower than searching for a "new astronomical object."

The workflow instead asks:

1. Is the observation statistically unusual?
2. Does the signal survive robustness checks?
3. Is there source-like evidence in the image?
4. Could artifacts or observing conditions explain it?
5. Could blending or neighboring sources explain it?
6. Does another infrared band reproduce the behavior?
7. Do external catalogs identify a known source?
8. Do known astronomical classes provide an explanation?
9. What remains unexplained?
10. What follow-up would most efficiently discriminate between explanations?

---

Evidence-First Pipeline

                    OBSERVED DATA
                         │
                         ▼
                   ┌───────────┐
                   │  UNUSUAL  │
                   └─────┬─────┘
                         │
                         ▼
              DATA-QUALITY VERIFIED
                         │
                         ▼
          ROBUSTNESS / OUTLIER CHECKS
                         │
                         ▼
        IMAGE / ARTIFACT / BLENDING CHECKS
                         │
                         ▼
        OBSERVATIONAL / SYSTEMATIC CHECKS
                         │
                         ▼
             CROSS-BAND COMPARISON
                         │
                         ▼
            CATALOG / KNOWN-CLASS CHECKS
                         │
                         ▼
               EVIDENCE SYNTHESIS
                         │
                         ▼
              FOLLOW-UP CANDIDATE
                         │
                         ▼
             OBSERVATIONAL CONFIRMATION

Each stage can weaken, strengthen, or leave unchanged the interpretation of a candidate.

Missing evidence is treated as unknown rather than automatically treated as negative evidence.

---

Dataset

NEOWISE

The primary data source is the NEOWISE Single-exposure Source Database, distributed through the NASA/IPAC Infrared Science Archive (IRSA).

NEOWISE provides repeated measurements in:

Band| Approximate wavelength
W1| 3.4 μm
W2| 4.6 μm

The data include information such as:

- sky position;
- photometry;
- photometric uncertainty;
- observation time;
- frame and quality information.

The public repository does not redistribute the full NEOWISE survey archive.

Researchers attempting reproduction should obtain the required survey data directly from the appropriate public archive.

---

Development Field

The development analysis used a small sky region approximately centered on:

RA  ≈ 179.9°–180.1°
Dec ≈ −0.1°–+0.1°

This field was used for methodological development and investigation.

«Important: The development field should not be interpreted as a survey-wide statistical sample.»

The results therefore do not establish a survey-wide prevalence of unusual sources or a formal survey-wide false-positive rate.

---

Analysis Snapshot

The analysis progressed through the following population reductions:

Stage| Population
Initial detections| 8,652
Quality/usability-filtered detections| 7,275
Source-like groups| 3,371
Sources with ≥10 observations| 175
Validator-positive groups| 6
Top unresolved Tier 1 candidates| 2

The two primary unresolved candidates are:

Source| Priority| Current interpretation
1978| Tier 1| Persistent unusual W1 behavior; highest follow-up priority
1974| Tier 1| Persistent unusual W1 behavior; high follow-up priority
233| Tier 2| Requires additional controls
1976| Tier 2| Requires additional controls
234| Control/sanity| Control-compatible
1979| Control/sanity| Control-compatible

No candidate is established as a new astronomical class, new physical phenomenon, or confirmed intrinsic astrophysical variability.

---

Candidate Evidence

Source 1978

Source 1978 showed persistent excess W1 variability relative to the uncertainty-aware baseline.

Key results:

Test| Result
Phase III W1 reduced χ²| ≈ 3.94
Phase VII W1 reduced χ²| ≈ 4.14
Leave-one-out minimum χ²| ≈ 2.80
W2 reduced χ²| ≈ 1.29
W1–W2 correlation| ≈ 0.303

The signal survived leave-one-observation-out testing.

Matched-control comparisons also found the source unusual within the tested comparison population.

A nearby NEOWISE source was identified, so blending or contamination remains a live hypothesis. Direct neighbor tracking did not strongly reproduce the target behavior, but this does not completely eliminate blending.

No valid SIMBAD identification was returned under the tested conditions. The Gaia query was unavailable.

These are unknown states, not evidence that the source has no counterpart.

Current interpretation

1978 is a high-priority unexplained candidate requiring additional investigation.

It is not a confirmed astrophysical discovery.

---

Source 1974

Source 1974 also showed persistent excess W1 variability.

Key results:

Test| Result
Phase III W1 reduced χ²| ≈ 3.03
Phase VII W1 reduced χ²| ≈ 3.05
Leave-one-out minimum χ²| ≈ 2.14
W2 reduced χ²| ≈ 1.17
W1–W2 correlation| ≈ 0.041

The W1 signal survived leave-one-observation-out testing.

Matched-control comparisons found the source unusual within the tested comparison population.

A nearby NEOWISE source was identified, so blending remains a possible explanation, although direct neighbor tracking did not strongly reproduce the target behavior.

No valid SIMBAD identification was returned under the tested conditions, and the Gaia query was unavailable.

Current interpretation

1974 is also a high-priority unexplained candidate requiring additional investigation.

It is not a confirmed astrophysical discovery.

---

What the Analysis Supports

The analysis supports the following statements:

- 1978 exhibits persistent unusual W1 temporal scatter under the operational analysis.
- 1974 exhibits persistent unusual W1 temporal scatter under the operational analysis.
- Both candidates retain their signal under leave-one-observation-out testing.
- Source-like image evidence exists at the expected positions.
- Tested neighboring-source behavior does not simply reproduce the candidate signals.
- Matched-control and observational/systematic checks do not provide a dominant explanation within the tested scope.
- The candidates therefore warrant further investigation.

---

What the Analysis Does Not Establish

The analysis does not establish:

- a new astronomical class;
- new physics;
- confirmed intrinsic astrophysical variability;
- absence of blending;
- absence of instrumental effects;
- absence of extraction/systematic effects;
- a formal survey-wide false-positive rate;
- periodicity;
- a complete external identity;
- observational confirmation.

The term "unexplained" has a deliberately limited meaning:

«The tested explanations did not adequately account for the observations available in this analysis.»

It does not mean:

«All possible explanations have been eliminated.»

---

Statistical Validation

The operational validator used:

Nobs ≥ 10
median W1 SNR proxy ≥ 5
W1 reduced χ² ≥ 2
good-quality fraction ≥ 0.90

The corrected validation produced:

6 positive groups
82 evaluable comparison groups

giving:

6 / 82 ≈ 7.32%

However, 7.32% is not interpreted as a false-positive rate.

There is no independently labeled non-variable ground-truth population in this analysis. The comparison population is therefore insufficient to support a formal false-positive estimate.

---

Methodological Architecture

The project deliberately separates three concepts:

1. Novelty

How unusual is the observed behavior relative to the analyzed population?

2. Evidence quality

How strong and independently supported is the evidence for the observed behavior?

3. Systematic risk

What observational, instrumental, extraction, blending, or catalog-related explanations remain plausible?

These concepts are not collapsed into a single arbitrary weighted score.

The broader workflow uses multiple evidence views:

- temporal behavior;
- photometric behavior;
- image evidence;
- observational/systematic evidence;
- external survey information;
- known-class similarity.

---

Machine Learning

Machine learning is used as a triage and diagnostic tool, not as a discovery oracle.

The project explores unsupervised novelty detection to help identify sources that deserve investigation.

ML output is therefore interpreted as:

statistical / feature-space unusualness

rather than:

astrophysical novelty

A source receiving a high anomaly score does not establish that the source is astrophysically unusual.

---

Public Notebook Collection

The public notebooks are a curated subset of the larger private research history.

01 — Initial NEOWISE Validation

"notebooks/public/01_initial_neowise_validation.ipynb"

Early validation of the NEOWISE candidate-selection and photometric workflow.

02 — Validation Population

"notebooks/public/02_validation_population.ipynb"

Construction and validation of the source population used for later comparisons.

03 — Image Analysis

"notebooks/public/03_image_analysis.ipynb"

Image-level investigation of candidate observations, source morphology, photometry, and artifact-related evidence.

04 — Phase IV ML Analysis

"notebooks/public/04_phase4_ml_analysis.ipynb"

Feature construction and exploratory unsupervised novelty analysis.

«ML is used for triage and diagnosis rather than confirmation of astrophysical novelty.»

05 — Analysis Freeze Audit

"notebooks/public/05_analysis_freeze_audit.ipynb"

Audit and provenance checks associated with the frozen analysis state.

---

Reusable Source Code

Reusable utilities are provided under:

src/
├── coordinates.py
├── statistics.py
└── validation.py

These modules provide:

- coordinate utilities;
- uncertainty-aware statistical calculations;
- variability-validation helpers extracted from the broader research workflow.

---

Quick Start

Clone the repository and install the required dependencies:

python -m pip install -r requirements.txt

Core dependencies include:

- NumPy
- Pandas
- SciPy
- Matplotlib
- Astropy
- PyArrow
- scikit-learn
- Requests

Large survey datasets are intentionally excluded from Git.

---

Expected Local Data Layout

A local reproduction environment may organize data as:

data/
├── raw/
├── processed/
└── results/

For example:

data/
└── raw/
    └── exp005_region_detections.csv

The exact input products required by individual historical experiments may differ.

«Do not commit raw NEOWISE survey files, large Parquet datasets, downloaded image archives, or the private research archive.»

---

Reproducibility Boundary

This repository is the public research-code and reproducibility layer.

The complete private research archive remains outside this repository.

The private archive contains:

- broader experimental history;
- intermediate artifacts;
- large data products;
- frozen evidence archives;
- working materials;
- other development files not required for the public code release.

Therefore, a fresh clone does not claim to reconstruct every historical experiment byte-for-byte.

The public release provides:

- reusable research code;
- selected research notebooks;
- dependency information;
- data-source documentation;
- methodology;
- provenance information;
- selected analysis context.

Some historical notebooks depend on intermediate data products generated during earlier experiments. Those products are not all included in this public repository.

Reproduction of a specific historical result may therefore require obtaining the corresponding public survey data and reconstructing the relevant intermediate products.

The frozen numerical conclusions documented in the research archive should not be silently regenerated or replaced merely to make the public repository appear more complete.

---

What Is Not Included

The following are intentionally excluded:

- full NEOWISE raw/survey datasets;
- large downloaded Parquet archives;
- downloaded image archives;
- private Google Drive paths;
- private research archives;
- personal working files;
- credentials or secrets;
- uncurated intermediate experiment outputs.

---

Analysis Freeze

The research analysis was frozen through EXP-056.

The freeze state recorded:

ANALYSIS_FREEZE_READY_WITH_DOCUMENTED_LIMITATIONS

The freeze included:

- candidate master;
- experiment provenance;
- freeze checks;
- synthesis outputs;
- source checksums;
- documented limitations and warnings.

The frozen state contains no hard failures and documented limitations are retained rather than hidden.

---

Scientific Limitations

Important limitations include:

- the development field is small and does not constitute a survey-wide population study;
- the validator-positive set is not a formal false-positive rate;
- matched-control comparisons are limited in size and field;
- the control population is not equivalent to a statistically complete survey population;
- missing external information is not treated as negative evidence;
- image-level source-like structure does not prove astrophysical variability;
- neighboring-source tests do not completely eliminate blending;
- unsupervised ML novelty is not equivalent to astrophysical novelty;
- exact ML ranks can be sensitive to population and feature sampling;
- temporal sampling limits interpretation of possible variability patterns;
- no candidate has observational confirmation in this analysis.

These limitations are part of the result rather than something hidden from the reader.

---

Evidence Language

The project uses conservative evidence labels.

Label| Meaning
CONFIRMED| Evidence supports confirmation under the defined analysis
SUPPORTED| Evidence provides substantive support
LIKELY| Evidence favors the interpretation but does not establish it
PLAUSIBLE| Consistent with the available evidence
INCONCLUSIVE| Evidence does not adequately distinguish between interpretations
REJECTED| Tested evidence does not support the interpretation
NOT TESTED| The relevant test was not performed
UNKNOWN| Required information is unavailable

These labels describe the strength or availability of evidence.

They should not be interpreted as probabilities unless an explicit statistical probability model has been constructed.

---

Data and Code Availability

The public repository contains the research code, selected notebooks, methodology documentation, and the canonical manuscript.

The underlying NEOWISE survey data and external catalog data are not redistributed.

Researchers should obtain those datasets from their respective public archives and reproduce only the data products required for the analysis.

---

Citation

If you use the research code or methodology, please cite this repository.

Citation metadata are provided in:

CITATION.cff

The accompanying manuscript is:

«Evidence-First Astronomical Anomaly Investigation»

---

License

This project is released under the MIT License.

See:

LICENSE

and

CITATION.cff

for licensing and citation metadata.

---

Research Principle

Detect → Verify → Explain

The purpose of anomaly detection is not to manufacture discoveries.

It is to identify observations that deserve investigation and then progressively test whether those observations survive increasingly strong alternative explanations.

«An anomaly is a hypothesis, not a discovery.»
