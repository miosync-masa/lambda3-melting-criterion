# Λ³-Based Lindemann Melting Criterion: A Unified Framework

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)

## Overview

This repository contains the complete computational framework and datasets supporting our manuscript on **physically-derived Lindemann melting criterion** based on the Λ³ (Lambda-cubed) energy density ratio theory.

We demonstrate that:
1. **Thermal softening** of elastic moduli follows a universal Λ³ exponential decay across 7 elemental metals
2. **Lindemann parameters** can be predicted from first principles using Debye-elastic theory with structure-dependent shear softening
3. **Nanoparticle melting** follows Gibbs-Thomson law (1/r) in the surface-dominated regime (r ≥ 3 nm)

**Manuscript**: [To be added upon acceptance]  
**Authors**: Masamichi Iizumi (Miosync Inc.) [Additional authors TBD]  
**Journal**: Nature Materials (under preparation)

---

## Repository Structure
```
lambda3-melting-criterion/
├── data/
│   ├── lindemann_parameters_literature.csv    # SI-3: Experimental δ_L with sources
│   ├── lindemann_physical_results.csv         # SI-4: Born-struct predictions
│   ├── material_summary.csv                   # SI-2: Λ³ fitted parameters
│   ├── nanoparticle_melting_gibbs_thomson.csv # SI-6: Size-dependent melting
│   └── thermal_properties.csv                 # SI-2: E(T) raw data
├── code/
│   ├── lindemann_predictor.py                 # SI-4: Lindemann criterion derivation
│   ├── gibbs_thomson_analysis.py              # SI-6: Nanoparticle analysis
│   └── lambda3_thermal_softening.py           # SI-2: Thermal softening fits
├── supplementary_information/
│   └── SI-1_Overview.md                       # SI-1: Methods and conventions
├── figures/                                    # Generated figures
├── README.md
├── requirements.txt
└── LICENSE
```

---

## Installation

### Requirements

- Python 3.10 or higher
- NumPy ≥ 1.21.0
- Pandas ≥ 1.3.0
- Matplotlib ≥ 3.4.0 (for visualization)
- SciPy ≥ 1.7.0 (for statistical analysis)

### Setup
```bash
# Clone the repository
git clone https://github.com/[username]/lambda3-melting-criterion.git
cd lambda3-melting-criterion

# Install dependencies
pip install -r requirements.txt
```

---

## Usage

### 1. Λ³ Thermal Softening Analysis (SI-2)

Fits Young's modulus temperature dependence E(T)/E₀ = exp[-λ_eff α ΔT] for 7 metals:
```bash
python code/lambda3_thermal_softening.py
```

**Output**:
- `Lambda3_individual_fits.png` - Material-by-material E(T) fits
- `Lambda3_fitted_parameters.csv` - Fitted λ_base and κ values

**Key Results**:
- All materials: Residual < 2.1%
- Mean accuracy: 1.0% (FCC), 1.3% (BCC), 1.2% (HCP)

---

### 2. Lindemann Criterion Prediction (SI-4)

Derives Lindemann parameters from Debye-elastic theory with Born-type shear softening:
```bash
python code/lindemann_predictor.py
```

**Output**:
- `lindemann_physical_results.csv` - Predicted vs experimental δ_L

**Key Results**:
- Baseline (harmonic): Mean error 69.9%
- Born-struct rule: Mean error **5.4%**

**Structure-dependent shear retention factors (f_G)**:
- BCC: 0.027
- FCC: 0.101
- HCP: 0.137

---

### 3. Nanoparticle Melting Analysis (SI-6)

Validates Gibbs-Thomson law (ΔT_m/T_m ∝ 1/r) for 7 metallic nanoparticles:
```bash
python code/gibbs_thomson_analysis.py
```

**Output**:
- `nanoparticle_piecewise_fit_english.png` - Size-dependent melting curves
- CSV with fitted parameters (a, b, R²)

**Key Results**:
- r ≥ 3 nm: Classical 1/r law (R² = 0.77–0.95)
- r < 3 nm: Edge effects require 1/r² correction

---

## Data Description

### `thermal_properties.csv`
Temperature-dependent material properties for 7 metals (Fe, W, Cu, Al, Ni, Ti, Mg):
- **Columns**: Material, T_K, T_C, alpha_1e-6, E_GPa, sigma_y_MPa, Notes
- **Source**: Compiled from ASM Handbook, ultrasonic measurements, DFT calculations

### `material_summary.csv`
Λ³ model parameters and material constants:
- **Columns**: Material, Structure, Tm_K, lambda_base, kappa, residual_percent, ...
- **Usage**: SI-2 Table 1

### `lindemann_parameters_literature.csv`
Experimental Lindemann melting parameters with literature sources:
- **Columns**: Material, Structure, Tm_K, delta_L_exp, delta_L_definition, Method, Source, ...
- **Usage**: SI-3 validation

### `lindemann_physical_results.csv`
Born-struct predictions vs experimental values:
- **Columns**: Material, Struct, δ_exp, fG_struct, δ_pred, Error_%
- **Usage**: SI-4 Table 1

### `nanoparticle_melting_gibbs_thomson.csv`
Size-dependent melting temperatures for nanoparticles:
- **Columns**: Material, Radius_nm, Tm_K, Tm_bulk_K, Tm_ratio, Delta_Tm_ratio, Source_Notes
- **Materials**: Au, Ag, Cu, Al, Fe, Pt, Pb
- **Usage**: SI-6 Gibbs-Thomson validation

---

## Key Findings

### 1. Thermal Softening (Λ³ Model)
All 7 metals follow exponential decay with material-specific parameters:
- **FCC metals**: Highest accuracy (0.7% mean residual)
- **Mg anomaly**: Extreme anharmonicity (κ ≈ 38) due to low Debye temperature

### 2. Lindemann Criterion
Structure-dependent shear softening at melting:
- **BCC**: 97% shear loss (f_G = 0.03)
- **FCC**: 90% shear loss (f_G = 0.09)
- **Prediction accuracy**: 5.4% mean error (vs 70% for harmonic approximation)

### 3. Nanoparticle Melting
- **Surface regime** (r ≥ 3 nm): Classical Gibbs-Thomson 1/r law validated
- **Edge regime** (r < 3 nm): Requires 1/r² correction
- **Physical basis**: Coordination-dependent Λ threshold

---

## Citation

If you use this code or data, please cite:
```bibtex
@article{iizumi2025lambda3,
  title={Physically-Derived Lindemann Melting Criterion via Λ³ Energy Density Ratio Theory},
  author={Iizumi, Masamichi and [TBD]},
  journal={Nature Materials},
  year={2025},
  note={In preparation}
}
```

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

Code and data are freely available for academic and commercial use with attribution.

---

## Contact

**Masamichi Iizumi**  
CEO, Miosync Inc.  
Email: [To be added]

**Issues & Questions**:  
Please open an issue on GitHub or contact the corresponding author.

---

## Acknowledgments

We thank our industrial partners (Nidec, NSK, JTEKT, NACHI-Fujikoshi) for experimental validation support and valuable discussions on practical applications.

Special thanks to the materials science community for providing high-quality experimental data used in this study.

---

## Version History

- **v1.0.0** (2025-XX-XX): Initial release supporting Nature Materials submission
  - Λ³ thermal softening model
  - Lindemann criterion derivation
  - Gibbs-Thomson nanoparticle analysis
  - Complete dataset compilation

---

**Last Updated**: 2025-11-14
