# Stellar Spectral Analysis Pipeline

[![Python](https://img.shields.io/badge/python-3.x-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![GitHub issues](https://img.shields.io/github/issues/yourusername/stellar-spectral-pipeline)](https://github.com/yourusername/stellar-spectral-pipeline/issues)
[![GitHub forks](https://img.shields.io/github/forks/yourusername/stellar-spectral-pipeline)](https://github.com/yourusername/stellar-spectral-pipeline/network)
[![GitHub stars](https://img.shields.io/github/stars/yourusername/stellar-spectral-pipeline)](https://github.com/yourusername/stellar-spectral-pipeline/stargazers)

---

## Overview

This repository contains a comprehensive pipeline for high-resolution stellar spectral analysis. It enables you to:

- Normalize observed spectra by continuum fitting  
- Determine spectral resolving power from absorption lines  
- Extract and filter synthetic spectra metadata  
- Degrade synthetic spectra to observed resolution  
- Apply rotational broadening to synthetic templates  
- Perform chi-square template matching for stellar parameter estimation  
- Convert MARCS atmosphere files to spectral synthesis compatible formats  
- Generate synthetic spectra automatically using Kurucz and MARCS model atmospheres with SPECTRUM  

This pipeline facilitates detailed stellar parameter determination and synthetic spectral fitting, crucial for astrophysics research on stars and exoplanet hosts.

---

## Table of Contents

- [Tasks](#tasks)  
  - [Task 1: Continuum Renormalization](#task-1-continuum-renormalization)  
  - [Task 2: Determining Resolving Power of Observed Spectrum](#task-2-determining-resolving-power-of-observed-spectrum)  
  - [Task 3: Synthetic Spectra Metadata Extraction and Filtering](#task-3-synthetic-spectra-metadata-extraction-and-filtering)  
  - [Task 4: Degrading Synthetic Spectra to Observed Resolving Power](#task-4-degrading-synthetic-spectra-to-observed-resolving-power)  
  - [Task 5: Rotational Broadening of Synthetic Spectra](#task-5-rotational-broadening-of-synthetic-spectra)  
  - [Task 6: Template Matching via Chi-Square Minimization](#task-6-template-matching-via-chi-square-minimization)  
  - [Automated Synthetic Spectra Generation Using Kurucz Supermodels](#automated-synthetic-spectra-generation-using-kurucz-supermodels)  
  - [MARCS Atmosphere Conversion](#marcs-atmosphere-conversion)  
  - [Synthetic Spectrum Generation Using SPECTRUM](#synthetic-spectrum-generation-using-spectrum)  
- [Requirements](#requirements)  
- [Usage](#usage)  
- [References](#references)  
- [Acknowledgments](#acknowledgments)  

---

## Tasks

### Task 1: Continuum Renormalization

- **Goal:** Normalize observed stellar spectra by estimating the continuum and removing absorption line effects.  
- **Method:**  
  - Perform spectral binning and sigma clipping to isolate continuum anchor points.  
  - Apply Savitzky-Golay filtering and spline fitting for continuum estimation.  
  - Normalize observed flux by the fitted continuum spline.  
  - Save continuum fits and residuals for diagnostics.  
- **Output:** Normalized observed spectrum and continuum spline.  
- **Code:** Python with `numpy`, `scipy`, `pandas`, and `matplotlib`.

---

### Task 2: Determining Resolving Power of Observed Spectrum

- **Goal:** Measure the resolving power \( R = \frac{\lambda}{\text{FWHM}} \) of the observed spectrum by fitting absorption lines.  
- **Method:**  
  - Segment wavelength range, detect absorption lines using negative peak detection.  
  - Fit Gaussian profiles to lines and compute FWHM.  
  - Calculate resolving power per line, then aggregate statistics.  
  - Visualize best Gaussian fits for quality control.  
- **Output:** CSV file with resolving power metrics and diagnostic plots.  
- **Tools:** Python with `scipy.optimize.curve_fit` and plotting.

---

### Task 3: Synthetic Spectra Metadata Extraction and Filtering

- **Goal:** Extract stellar parameters from synthetic spectrum filenames and select spectra within desired parameter ranges.  
- **Method:**  
  - Parse filenames for [M/H], \( T_{\mathrm{eff}} \), \(\log g\), and \(v \sin i\).  
  - Store metadata in a CSV catalog.  
  - Filter based on user-supplied parameter ranges.  
  - Normalize synthetic spectra using continuum spline from Task 1.  
  - Compute synthetic spectra resolving power as in Task 2.  
- **Output:** Filtered synthetic spectra metadata CSV and normalized flux arrays.

---

### Task 4: Degrading Synthetic Spectra to Observed Resolving Power

- **Goal:** Match synthetic spectra resolution to that of the observed spectrum.  
- **Method:**  
  - Calculate required Gaussian convolution kernel width to degrade synthetic resolution to observed.  
  - Convolve synthetic flux arrays with Gaussian kernel.  
  - Normalize degraded spectra with percentile-based continuum fitting.  
  - Interpolate degraded spectra onto observed wavelength grid.  
  - Save degraded spectra and metadata.  
  - Visualize degraded vs observed spectra.  
- **Tools:** `scipy.signal.fftconvolve`, `scipy.interpolate`, `numpy`.

---

### Task 5: Rotational Broadening of Synthetic Spectra

- **Goal:** Apply rotational broadening to synthetic spectra to simulate stellar rotation effects.  
- **Method:**  
  - Use limb-darkened rotational broadening kernel with limb-darkening coefficient \(\epsilon = 0.6\).  
  - Compute velocity grid and perform convolution in velocity space using `fftconvolve`.  
  - For \(v \sin i = 0\), return unaltered spectra.  
  - Store rotationally broadened templates for matching.  
  - Visualize templates vs observed spectra.  
- **Tools:** `scipy.signal`, `numpy`.

---

### Task 6: Template Matching via Chi-Square Minimization

- **Goal:** Determine best-fit synthetic spectral templates by minimizing chi-square differences with observed spectra.  
- **Method:**  
  - Compute \(\chi^2 = \sum (F_{\mathrm{obs}} - F_{\mathrm{syn}})^2\) for each template.  
  - Skip spectra with flux length mismatches.  
  - Rank and save chi-square values.  
  - Plot observed spectrum overlaid with best and second-best matches.  
- **Output:** CSV of chi-square results and matching plots.

---

### Automated Synthetic Spectra Generation Using Kurucz Supermodels

Kurucz supermodels are stellar atmosphere grids widely used for spectral synthesis.

- **Goal:** Automate the generation of synthetic spectra using Kurucz supermodels and SPECTRUM.  
- **Features:**  
  - Grid generation over realistic ranges of \( T_{\mathrm{eff}} \), \(\log g\), [M/H], and \(v \sin i\).  
  - Empirical microturbulence velocity estimation based on \( T_{\mathrm{eff}} \) and \(\log g\).  
  - Extraction of atmospheric models from Kurucz supermodel files via `selectmod`.  
  - Synthetic spectra computed with SPECTRUM for user-defined wavelength range and resolution.  
  - Rotational broadening applied via `avsini`.  
  - Robust checking of required binaries (`selectmod`, `spectrum`, `avsini`) and inputs (`luke.lst`, `Supermodels/`).  
  - Organized output directories for extracted models and spectra.  
  - Error handling and user confirmation prompts.  

---

### MARCS Atmosphere Conversion

MARCS models are 1D LTE hydrostatic atmospheres for F, G, K, M stars.

- **Goal:** Convert MARCS `.krz` atmospheric structure files into Kurucz-compatible `.mod` files.  
- **Steps:**  
  - Parse MARCS `.krz` files for stellar parameters and layer data.  
  - Compute gas pressure from electron density and temperature using ideal gas law.  
  - Include microturbulence velocity from filename or default.  
  - Output `.mod` files with structure needed by SPECTRUM.  
  - Batch conversion over all `.krz` files in a directory.  
- **Output:** `.mod` files saved in a user-specified directory.  

---

### Synthetic Spectrum Generation Using SPECTRUM

- **Goal:** Automate synthetic spectrum generation from `.mod` atmosphere files using SPECTRUM.  
- **Features:**  
  - Copy `.mod` files into SPECTRUM executable directory.  
  - Parse stellar parameters from `.mod` filenames.  
  - Validate parameters to ensure safe operation of SPECTRUM.  
  - Run SPECTRUM in batch mode for a given wavelength range and resolution.  
  - Collect output `.spc` files into a dedicated output folder.  
  - Handle timeouts and runtime errors gracefully.  
  - Clean up intermediate files post-processing.  

---

## Requirements

- Python 3.x  
- Python Packages:  
  - `numpy`  
  - `scipy`  
  - `pandas`  
  - `matplotlib`  
  - `astropy`  
  - `tqdm` (for progress bars)  
- SPECTRUM software binaries (`selectmod`, `spectrum`, `avsini`) compiled and executable  
- Kurucz Supermodel files in `Supermodels/` directory  
- MARCS `.krz` atmosphere model files (if using MARCS tasks)  
- Line list file `luke.lst` compatible with SPECTRUM  

---

## Usage

### Running Tasks 1 to 6 (Spectral analysis and template matching)

```bash
# Run continuum normalization script (Task 1)
python continuum_normalization.py

# Run resolving power determination (Task 2)
python resolve_power.py

# Extract and filter synthetic spectra metadata (Task 3)
python synthetic_metadata.py

# Degrade synthetic spectra to observed resolution (Task 4)
python degrade_spectra.py

# Apply rotational broadening to synthetic spectra (Task 5)
python rotational_broadening.py

# Perform chi-square template matching (Task 6)
python chi_square_matching.py

