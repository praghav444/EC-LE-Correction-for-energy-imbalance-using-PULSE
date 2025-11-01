# PULSE: A Method to Correct Latent Heat Flux Measurements for Energy Imbalance at Eddy Covariance Sites

**Authors:**  
Pushpendra Raghav¹, Mukesh Kumar¹  
¹ Department of Civil, Construction, and Environmental Engineering, The University of Alabama  

---

## 📘 Overview

**PULSE** (Potential Underlying Water Use Efficiency-based Latent Heat Correction) provides a novel approach to correct **latent heat (LE) fluxes** for **surface energy imbalance** observed at **eddy covariance (EC)** flux tower sites.

This repository contains an R-based tutorial demonstrating the application of the PULSE method on a sample FLUXNET2015-style dataset. The workflow estimates corrected LE fluxes using the concept of **potential underlying water-use efficiency (uWUEₚ)** — an ecohydrological benchmark that helps identify and adjust periods of likely LE underestimation.

---

## 📂 Repository Contents

| File | Description |
|------|--------------|
| `run_PULSE_tutorial.R` | Main R script demonstrating PULSE implementation step-by-step |
| `sample_data_FR-Pue_2000_2014.csv` | Example dataset (subset of FLUXNET site FR-Pue) |
| `README_PULSE_Tutorial.md` | This instruction and documentation file |
| `LICENSE.txt` | License terms for reuse and citation |

---

## 🧭 Getting Started

### 1. Requirements

- **R version:** ≥ 4.3.0  
- Tested on Windows, macOS, and Linux environments.

### 2. Required Packages

The script automatically installs any missing dependencies.  
Alternatively, install them manually using:

```r
install.packages(c("lubridate", "dplyr", "data.table", "quantreg", "bigleaf", "signal"))
