# fMRIPrep to CONN Batch Processing with QC and Scrubbing

A MATLAB script for loading fMRIPrep-preprocessed fMRI data into the CONN toolbox with automated quality control (QC) and motion scrubbing.

## Quick Start

### 1. Clone the Repository
```bash
git clone <repository-url>
cd <repository-name>
```

### 2. Run the MATLAB Script
```bash
# Open MATLAB
matlab

# In MATLAB GUI:
# Navigate to the cloned repository directory using the file browser
# Then run the script:
batch_fmriprep_conn
```

### 3. Load the Project in CONN
```bash
# After the script completes successfully, launch CONN from MATLAB terminal:
conn

# In the CONN GUI:
# File > Open > Navigate to your project_dir
# Select: conn_project.mat
```

## Overview

This script automates the process of:
- Loading fMRIPrep derivatives into CONN toolbox
- Extracting and organizing anatomical ROIs (Grey Matter, White Matter, CSF)
- Processing motion confounds and creating scrubbing regressors
- Computing quality control metrics (FD, DVARS)
- Setting up denoising parameters for functional connectivity analysis

## Features

- ✅ **Multi-subject processing** - Batch process multiple subjects from CSV file
- ✅ **Multi-session support** - Handles single or multiple scanning sessions per subject
- ✅ **Multi-run concatenation** - Properly concatenates runs within sessions
- ✅ **Motion scrubbing** - Creates spike regressors for high-motion timepoints
- ✅ **Quality control** - Computes and saves FD, DVARS, and scrubbing statistics
- ✅ **fMRIPrep integration** - Uses fMRIPrep's motion_outlier regressors when available
- ✅ **Flexible thresholds** - Customizable FD and DVARS thresholds

## Requirements

### Software
- MATLAB (R2018b or later recommended)
- [CONN Toolbox](https://web.conn-toolbox.org/) (v20.b or later)
- [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) (required by CONN)

### Data Structure
Your data should be organized in BIDS-compliant fMRIPrep output format:

**Single-session data:**
```
fmriprep_dir/
├── sub-XXXX/
│   ├── anat/
│   │   ├── sub-XXXX_space-MNI152NLin2009cAsym_desc-preproc_T1w.nii.gz
│   │   ├── sub-XXXX_space-MNI152NLin2009cAsym_label-GM_probseg.nii.gz
│   │   ├── sub-XXXX_space-MNI152NLin2009cAsym_label-WM_probseg.nii.gz
│   │   └── sub-XXXX_space-MNI152NLin2009cAsym_label-CSF_probseg.nii.gz
│   └── func/
│       ├── sub-XXXX_task-XXX_run-01_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz
│       ├── sub-XXXX_task-XXX_run-01_desc-confounds_timeseries.tsv
│       ├── sub-XXXX_task-XXX_run-02_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz
│       └── sub-XXXX_task-XXX_run-02_desc-confounds_timeseries.tsv
└── sub-YYYY/
    └── ...
```

**Multi-session data (BIDS format):**
```
fmriprep_dir/
├── sub-XXXX/
│   ├── anat/
│   │   ├── sub-XXXX_space-MNI152NLin2009cAsym_desc-preproc_T1w.nii.gz
│   │   ├── sub-XXXX_space-MNI152NLin2009cAsym_label-GM_probseg.nii.gz
│   │   ├── sub-XXXX_space-MNI152NLin2009cAsym_label-WM_probseg.nii.gz
│   │   └── sub-XXXX_space-MNI152NLin2009cAsym_label-CSF_probseg.nii.gz
│   ├── ses-01/
│   │   └── func/
│   │       ├── sub-XXXX_ses-01_task-rest_run-01_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz
│   │       ├── sub-XXXX_ses-01_task-rest_run-01_desc-confounds_timeseries.tsv
│   │       ├── sub-XXXX_ses-01_task-rest_run-02_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz
│   │       └── sub-XXXX_ses-01_task-rest_run-02_desc-confounds_timeseries.tsv
│   └── ses-02/
│       └── func/
│           ├── sub-XXXX_ses-02_task-rest_run-01_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz
│           └── sub-XXXX_ses-02_task-rest_run-01_desc-confounds_timeseries.tsv
└── sub-YYYY/
    └── ...
```

The script automatically detects whether your data has sessions and processes accordingly.

## Setup

### 1. Configure Paths

Edit the following variables in `batch_fmriprep_conn.m`:

```matlab
fmriprep_dir = '/path/to/your/fmriprep/derivatives';
project_dir  = '/path/to/your/conn/project';
participant_file = fullfile(project_dir, 'participant_ids.csv');
TR = 2.0;  % Repetition time in seconds
```

### 2. Create Participant List

Create a `participant_ids.csv` file in your project directory:

```csv
participant_id
sub-1656
sub-1700
sub-1708
```

An example file is included in this repository.

### 3. Adjust QC Thresholds (Optional)

Modify quality control parameters as needed:

```matlab
FD_threshold = 0.5;      % Framewise displacement threshold (mm)
DVARS_threshold = 1.5;   % DVARS threshold (standardized)
```

**Common thresholds:**
- **Conservative:** FD < 0.2 mm
- **Moderate:** FD < 0.5 mm (default)
- **Liberal:** FD < 0.9 mm

## Usage

### Running the Script

1. Open MATLAB
2. Navigate to the script directory
3. Run the script:

```matlab
batch_fmriprep_conn
```

### Expected Output

The script will create:

```
project_dir/
├── conn_project.mat           # Main CONN project file
├── batch_prerun.mat           # Saved batch structure (for debugging)
├── qc_summary.mat             # Quality control metrics
└── (various CONN output files and folders)
```

For each session, motion and scrubbing files will be saved in the functional directory:
```
fmriprep_dir/sub-XXXX/ses-XX/func/
├── rp_sub-XXXX_session01.txt     # Motion parameters (6 DOF)
└── scrub_sub-XXXX_session01.txt  # Scrubbing regressors (spike regressors)
```

### Console Output

The script provides detailed progress information:

```
=== STARTING DATA COLLECTION ===

--- Subject: sub-1656 ---
T1w: sub-1656_space-MNI152NLin2009cAsym_desc-preproc_T1w.nii.gz
Number of sessions: 1

  Session 1: ses-01
  Found 3 functional files
    Run 1: sub-1656_ses-01_task-rest_run-01_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz
      Found confounds: sub-1656_ses-01_task-rest_run-01_desc-confounds_timeseries.tsv
      QC: FD=0.234 mm, DVARS=1.123, Scrubbed=12 (4.3%)
    ...
  Created motion file: rp_sub-1656_session01.txt (840 timepoints)
  Created scrubbing file: scrub_sub-1656_session01.txt (840 timepoints, 12 regressors)

=== QC SUMMARY ===
Mean FD across all runs: 0.287 ± 0.156 mm
Mean DVARS across all runs: 1.245 ± 0.234
Mean % scrubbed: 5.2% ± 3.1%

✓ SUCCESS! Project created: /path/to/project/conn_project.mat
✓ QC metrics saved to: /path/to/project/qc_summary.mat
```

## Quality Control Metrics

### Framewise Displacement (FD)
- Measures head motion between consecutive volumes
- Sum of absolute displacements (translation + rotation × 50mm)
- **Interpretation:** FD > 0.5 mm indicates significant motion

### DVARS
- Measures signal change between consecutive volumes
- Standardized root mean square change in BOLD signal
- **Interpretation:** High DVARS indicates artifacts or rapid signal changes

### Scrubbing
- High-motion timepoints are flagged for regression
- Each flagged timepoint gets a separate regressor (spike regressor)
- **Recommendation:** Exclude sessions with >20% scrubbed volumes

## Denoising Strategy

The script configures CONN with the following denoising parameters:

- **Bandpass filter:** 0.01 - 0.1 Hz (default for resting-state)
- **Detrending:** Linear
- **Despiking:** Off (already performed by fMRIPrep)
- **Confound regression:**
  - 6 motion parameters (translation + rotation)
  - Scrubbing regressors (spike regressors for outlier timepoints)
  - White matter signal
  - CSF signal

## Troubleshooting

### Common Issues

**1. "participant_ids.csv not found"**
- Ensure the CSV file exists in your project directory
- Check the `participant_file` path variable

**2. "T1w not found for sub-XXXX"**
- Verify anatomical files exist in `sub-XXXX/anat/`
- Check file naming matches fMRIPrep conventions

**3. "Index exceeds the number of array elements"**
- This was a bug in earlier versions - ensure you have the fixed version
- Check that covariate initialization uses `cell(2, 1)` not `cell(2, nsub)`

**4. "No functional files found"**
- Verify functional data exists in expected session/func directories
- Check that files match the space (MNI152NLin2009cAsym) and format (.nii.gz)

**5. High percentage of scrubbed volumes (>20%)**
- Consider excluding the subject/session from analysis
- Review original data quality
- Adjust FD/DVARS thresholds if appropriate for your study

### Debug Mode

If the script fails, it saves debugging information:

```matlab
load(fullfile(project_dir, 'debug_batch_error.mat'));
```

This contains the `batch` structure and error details (`ME`).

## Script Components

### Main Functions

1. **Setup Anatomical ROIs** - Loads tissue probability maps from fMRIPrep
2. **Main Processing Loop** - Iterates through subjects, sessions, and runs
3. **Confound Extraction** - Extracts motion parameters and creates scrubbing regressors
4. **Concatenation** - Properly aligns multi-run data in time
5. **QC Computation** - Calculates and saves quality metrics
6. **CONN Batch Setup** - Configures CONN project structure
7. **Denoising Configuration** - Sets up filtering and confound regression

### Helper Functions

- `concatenate_scrubbing_regressors()` - Aligns scrubbing regressors across runs
- `extract_confounds_with_scrubbing()` - Extracts motion and QC metrics from TSV
- `convert_to_numeric()` - Handles mixed data types in confounds TSV
- `compute_FD()` - Calculates framewise displacement from motion parameters

## Citation

If you use this script, please cite:

- **CONN Toolbox:** Whitfield-Gabrieli S, Nieto-Castanon A. Conn: A Functional Connectivity Toolbox for Correlated and Anticorrelated Brain Networks. Brain Connect. 2012;2(3):125-141.

- **fMRIPrep:** Esteban O, Markiewicz CJ, Blair RW, et al. fMRIPrep: a robust preprocessing pipeline for functional MRI. Nat Methods. 2019;16(1):111-116.

- **Motion scrubbing:** Power JD, Barnes KA, Snyder AZ, Schlaggar BL, Petersen SE. Spurious but systematic correlations in functional connectivity MRI networks arise from subject motion. Neuroimage. 2012;59(3):2142-2154.

## References

- [CONN Toolbox Documentation](https://web.conn-toolbox.org/resources/documentation)
- [fMRIPrep Documentation](https://fmriprep.org/)
- [Power et al. (2012) - FD computation](https://doi.org/10.1016/j.neuroimage.2011.10.018)
- [Ciric et al. (2017) - Benchmarking denoising strategies](https://doi.org/10.1016/j.neuroimage.2017.03.020)

## License

This script is provided as-is for research purposes. Please ensure compliance with licenses for CONN Toolbox and SPM12.

## Author & Support

For questions or issues, please open an issue on the GitHub repository.

---

**Last Updated:** January 2026
