# Stellar Spectral Analysis Pipeline

[![Python](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![GitHub issues](https://img.shields.io/github/issues/yourusername/stellar-spectral-pipeline)](https://github.com/yourusername/stellar-spectral-pipeline/issues)
[![GitHub forks](https://img.shields.io/github/forks/yourusername/stellar-spectral-pipeline)](https://github.com/yourusername/stellar-spectral-pipeline/network)
[![GitHub stars](https://img.shields.io/github/stars/yourusername/stellar-spectral-pipeline)](https://github.com/yourusername/stellar-spectral-pipeline/stargazers)

---

## Overview

This repository hosts a robust pipeline for high-resolution stellar spectral analysis, designed to support detailed stellar parameter determination and synthetic spectral fitting—key tasks in astrophysics research, particularly for studying stars and exoplanet host systems. The pipeline integrates multiple scientifically validated tasks, leveraging both observed and synthetic spectra to derive precise stellar properties. It includes the following capabilities:

1. **Continuum Renormalization**: Normalizes observed spectra by fitting and removing the continuum, ensuring consistency across instruments and enabling accurate parameter estimation.
2. **Determining Resolving Power of Observed Spectrum**: Quantifies the spectral resolving power \( R = \frac{\lambda}{\text{FWHM}} \) using absorption line fits, critical for instrument characterization and synthetic spectra comparison.
3. **Synthetic Spectra Metadata Extraction and Filtering**: Parses and filters synthetic spectra metadata, streamlining template selection for analysis.
4. **Degrading Synthetic Spectra to Observed Resolving Power**: Adjusts synthetic spectra resolution to match observed spectra, ensuring fair comparisons.
5. **Rotational Broadening of Synthetic Spectra**: Applies rotational broadening to account for stellar rotation effects on line profiles.
6. **Template Matching via Chi-Square Minimization**: Identifies best-fit synthetic templates through chi-square minimization, a standard method for stellar parameter estimation.
7. **Automated Synthetic Spectra Generation Using Kurucz Supermodels**: Generates synthetic spectra automatically using Kurucz models and SPECTRUM, providing a robust template library.
8. **MARCS Atmosphere Conversion**: Converts MARCS atmospheric models into a format compatible with spectral synthesis tools.
9. **Synthetic Spectrum Generation Using SPECTRUM**: Automates spectrum synthesis from atmospheric models in batch mode.

This pipeline is built to handle high-precision astrophysical data and is suitable for production-level research applications.

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

- **Goal**: Normalize observed stellar spectra by estimating and removing the continuum to facilitate cross-instrument comparisons and parameter estimation.
- **Method**:
  - Bin the spectrum and apply sigma clipping to isolate continuum anchor points.
  - Smooth the spectrum using a Savitzky-Golay filter (window size: 51, polynomial order: 3).
  - Fit a cubic spline to the filtered continuum points.
  - Divide the observed flux by the fitted spline to normalize the spectrum.
  - Save continuum fits and residuals for diagnostic analysis.
- **Output**: Normalized spectrum and continuum spline files.
- **Tools**: Python with `numpy`, `scipy`, `pandas`, `matplotlib`.

---

### Task 2: Determining Resolving Power of Observed Spectrum

- **Goal**: Measure the resolving power \( R = \frac{\lambda}{\text{FWHM}} \) of observed spectra by analyzing absorption lines.
- **Method**:
  - Segment the wavelength range into intervals for line detection.
  - Identify absorption lines using a negative peak detection algorithm with a minimum depth threshold.
  - Fit Gaussian profiles to detected lines via least-squares optimization.
  - Calculate FWHM from Gaussian parameters and derive \( R \) per line.
  - Compute aggregate statistics (mean, standard deviation) and visualize fits.
- **Output**: CSV file with resolving power metrics and diagnostic plots.
- **Tools**: Python with `scipy.optimize.curve_fit`, `matplotlib`.

---

### Task 3: Synthetic Spectra Metadata Extraction and Filtering

- **Goal**: Extract stellar parameters from synthetic spectra filenames and filter based on user-defined ranges.
- **Method**:
  - Use regular expressions to parse filenames for [M/H], \( T_{\mathrm{eff}} \), \(\log g\), and \(v \sin i\).
  - Store metadata in a CSV catalog.
  - Apply filters based on user-specified parameter ranges.
  - Normalize synthetic spectra using the continuum spline from Task 1.
  - Compute resolving power as per Task 2.
- **Output**: Filtered metadata CSV and normalized flux arrays.
- **Tools**: Python with `pandas`, `numpy`, `re`.

---

### Task 4: Degrading Synthetic Spectra to Observed Resolving Power

- **Goal**: Adjust synthetic spectra resolution to match that of observed spectra.
- **Method**:
  - Calculate the Gaussian convolution kernel width required to degrade synthetic resolution.
  - Convolve synthetic flux arrays with the kernel using Fast Fourier Transform (FFT).
  - Normalize degraded spectra via percentile-based continuum fitting.
  - Interpolate onto the observed wavelength grid using cubic splines.
  - Save degraded spectra and metadata; visualize against observed spectra.
- **Tools**: `scipy.signal.fftconvolve`, `scipy.interpolate`, `numpy`, `matplotlib`.

---

### Task 5: Rotational Broadening of Synthetic Spectra

- **Goal**: Simulate stellar rotation effects by applying rotational broadening to synthetic spectra.
- **Method**:
  - Employ a limb-darkened rotational broadening kernel (\(\epsilon = 0.6\)).
  - Compute a velocity grid based on wavelength range and maximum \(v \sin i\).
  - Convolve in velocity space using `fftconvolve`.
  - Return unaltered spectra if \(v \sin i = 0\).
  - Store broadened templates and visualize against observed spectra.
- **Tools**: `scipy.signal`, `numpy`, `matplotlib`.

---

### Task 6: Template Matching via Chi-Square Minimization

- **Goal**: Identify best-fit synthetic templates by minimizing chi-square differences with observed spectra.
- **Method**:
  - Calculate \(\chi^2 = \sum (F_{\mathrm{obs}} - F_{\mathrm{syn}})^2\) for each template.
  - Exclude templates with mismatched flux lengths.
  - Rank templates by \(\chi^2\), save results, and plot the top two matches overlaid with the observed spectrum.
- **Output**: CSV of \(\chi^2\) results and matching plots.
- **Tools**: Python with `numpy`, `pandas`, `matplotlib`.

---

### Automated Synthetic Spectra Generation Using Kurucz Supermodels

Kurucz supermodels provide a widely-used grid of stellar atmospheres for spectral synthesis.

- **Goal**: Automate the creation of synthetic spectra using Kurucz supermodels and SPECTRUM.
- **Features**:
  - Generate grids over realistic ranges of \( T_{\mathrm{eff}} \), \(\log g\), [M/H], and \(v \sin i\).
  - Estimate microturbulence velocity empirically from \( T_{\mathrm{eff}} \) and \(\log g\).
  - Extract models using the `selectmod` binary.
  - Compute spectra with SPECTRUM for specified wavelength ranges and resolutions.
  - Apply rotational broadening with `avsini`.
  - Verify presence of binaries (`selectmod`, `spectrum`, `avsini`) and inputs (`luke.lst`, `Supermodels/`).
  - Organize outputs into directories; include error handling and user prompts.

---

### MARCS Atmosphere Conversion

MARCS models are 1D LTE hydrostatic atmospheres tailored for F, G, K, and M stars.

- **Goal**: Convert MARCS `.krz` files into Kurucz-compatible `.mod` files for spectral synthesis.
- **Method**:
  - Parse `.krz` files to extract stellar parameters and atmospheric layer data.
  - Compute gas pressure using the ideal gas law from electron density and temperature.
  - Incorporate microturbulence velocity from filenames or a default value.
  - Generate `.mod` files compatible with SPECTRUM.
  - Perform batch conversion across a directory of `.krz` files.
- **Output**: `.mod` files in a user-specified directory.
- **Tools**: Python with `numpy`, `pandas`.

---

### Synthetic Spectrum Generation Using SPECTRUM

- **Goal**: Automate synthetic spectrum generation from `.mod` files using SPECTRUM.
- **Features**:
  - Transfer `.mod` files to the SPECTRUM executable directory.
  - Parse and validate stellar parameters from filenames.
  - Execute SPECTRUM in batch mode for specified wavelength ranges and resolutions.
  - Collect `.spc` output files in a dedicated folder.
  - Handle runtime errors and timeouts; clean up intermediate files.

---

## Requirements

- **Python**: 3.8 or later
- **Python Packages**:
  - `numpy>=1.20`
  - `scipy>=1.7`
  - `pandas>=1.3`
  - `matplotlib>=3.4`
  - `astropy>=5.0`
  - `tqdm>=4.60` (progress bars)
- **External Software**:
  - SPECTRUM binaries (`selectmod`, `spectrum`, `avsini`), compiled and executable
- **Data Files**:
  - Kurucz Supermodel files in `Supermodels/` directory
  - MARCS `.krz` files (for MARCS tasks)
  - SPECTRUM-compatible line list `luke.lst`

---

## Usage

### Running Tasks 1 to 6 (Spectral Analysis and Template Matching)

Each task is executed via a dedicated Python script. Configure input/output paths in each script prior to execution.

```bash
# Task 1: Continuum normalization
python continuum_normalization.py

# Task 2: Resolving power determination
python resolve_power.py

# Task 3: Synthetic metadata extraction and filtering
python synthetic_metadata.py

# Task 4: Degrade synthetic spectra
python degrade_spectra.py

# Task 5: Apply rotational broadening
python rotational_broadening.py

# Task 6: Chi-square template matching
python chi_square_matching.py

# Generate synthetic spectra with Kurucz supermodels
python generate_kurucz_spectra.py

# Convert MARCS atmospheres to .mod files
python convert_marcs_to_mod.py

# Generate spectra from .mod files using SPECTRUM
python generate_spectrum_from_mod.py ```

References

Kurucz, R. L. (1993). ATLAS9 Stellar Atmosphere Programs and 2 km/s Grid. CD-ROM No. 13.
Gustafsson, B., et al. (2008). A grid of MARCS model atmospheres for late-type stars. Astronomy & Astrophysics, 486(3), 951-970.
Sbordone, L., et al. (2004). The SPECTRUM program for the synthesis of stellar spectra. Memorie della Societa Astronomica Italiana, 75, 416.

