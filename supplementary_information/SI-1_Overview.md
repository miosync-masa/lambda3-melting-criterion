# SUPPLEMENTARY INFORMATION SECTION SI-1

# Overview, Methodology, and Conventions

---

## SI-1.1 Organization of Supplementary Information

This Supplementary Information (SI) document provides comprehensive technical details, data sources, fitting procedures, and statistical validation supporting the main text. The SI is organized as follows:

**SI-2: Thermal Softening (Λ³ Model) – Data Sources and Fitting Procedures**
- Complete temperature-dependent Young's modulus E(T) data for seven metals (Fe, W, Cu, Al, Ni, Ti, Mg)
- Detailed fitting methodology for Λ³ parameters (λ_base, κ)
- Material-by-material validation with residual analysis
- *Correspondence*: Main text Section 3, Figure 2, Table 1

**SI-3: Lindemann Parameters – Literature Compilation and Normalization**
- Experimental and theoretical Lindemann melting parameters (δ_L) for seven metals
- Data provenance with complete literature references
- Normalization procedures (lattice constant basis → nearest-neighbor basis)
- Structure-dependent trends (BCC, FCC, HCP)
- *Correspondence*: Main text Section 4, Table 3

**SI-4: Born-Type Shear Collapse and Lindemann Predictions**
- Derivation of predicted versus required shear retention factors (f_G^pred, f_G^req)
- Structure-dependent clustering analysis
- Born-struct rule validation (mean absolute error 5.4%)
- Temperature evaluation sensitivity
- *Correspondence*: Main text Section 4, Table 3

**SI-5: EDR Regression – Statistical Diagnostics and Validation**
- Energy-density ratio (EDR) calculation at representative temperatures
- Through-origin regression analysis (δ_L = s√EDR)
- Leave-one-out cross-validation (LOO-CV)
- Bootstrap confidence intervals and permutation tests
- Model selection criteria and residual diagnostics
- *Correspondence*: Main text Section 5, Figure 4

**SI-6: Nanoparticle Melting-Point Depression – Complete Dataset and Analysis**
- Experimental melting data for seven metallic nanoparticles (Au, Ag, Cu, Al, Fe, Pt, Pb)
- Gibbs–Thomson relation validation (1/r law for r ≥ 3 nm)
- Size-dependent coordination reduction and EDR amplification
- Piecewise fitting procedure and statistical validation
- *Correspondence*: Main text Section 6, Figure 5

---

## SI-1.2 Computational Methods and Software Environment

All numerical calculations, nonlinear fitting, and statistical analyses were performed using Python with the following environment:

**Software versions:**
- Python: 3.10 or higher (tested on 3.10–3.12)
- NumPy: ≥1.21.0 (array operations, mathematical functions)
- Pandas: ≥1.3.0 (data manipulation, CSV I/O)
- Matplotlib: ≥3.4.0 (visualization)
- SciPy: ≥1.7.0 (optimization, statistical analysis)

**Core computational methods:**

*Nonlinear least-squares fitting (SI-2):*
Λ³ thermal-softening parameters (λ_base, κ) were extracted via the Levenberg–Marquardt algorithm implemented in `scipy.optimize.curve_fit`. Parameter bounds were imposed to ensure physical validity (λ_base ∈ [1, 200], κ ∈ [0.01, 50]). Convergence criterion: relative tolerance 10⁻⁸ in parameter space.

*Binary search inversion (SI-4):*
Required shear retention factors f_G^req were determined by numerically inverting the Debye–elastic Lindemann formula. A binary search on the interval f_G ∈ [10⁻⁶, 1.0] was iterated until |δ_calc − δ_exp| < 10⁻⁶. Typical convergence required 15–20 iterations per material.

*Physical constants:*
All calculations used SI units with the following fundamental constants:
- Boltzmann constant: k_B = 1.380649×10⁻²³ J/K
- Atomic mass unit: u = 1.66053906660×10⁻²⁷ kg
- Debye approximation for phonon density of states with isotropic elasticity

*Numerical stability:*
For materials with extreme thermal softening (Al, Mg), exponential overflow in the Λ³ model was mitigated by clipping the argument of exp[−λ_eff α ΔT] to a minimum of −40.

**Reproducibility:**
Complete Python scripts for all calculations, including data processing, fitting, inversion, and figure generation, are provided in the code repository (Section SI-1.6). All random-number-dependent procedures use fixed seeds for reproducibility.

---

## SI-1.3 Statistical Validation Framework

Multiple complementary statistical methods were employed to validate model predictions and assess uncertainty:

**Leave-one-out cross-validation (LOO-CV) (Sections SI-2, SI-5, SI-6):**
Predictive stability was assessed by removing each material (or data point) in turn, refitting the model on the remaining N−1 samples, and computing out-of-sample prediction error. LOO root-mean-square error (RMSE_LOO) quantifies generalization performance.

**Permutation tests (Sections SI-3, SI-5):**
Assumption-free significance testing was performed via exact permutation of material labels. Two-sided p-values quantify the probability of observing the test statistic under the null hypothesis without parametric distributional assumptions.

**Model selection criteria (SI-5):**
Through-origin versus intercept models were compared via Akaike Information Criterion (AIC) and uncentered R². Physical theory (δ ∝ √EDR with zero intercept) guided model selection, supplemented by statistical diagnostics.

**Residual diagnostics (SI-5):**
Normality of residuals was assessed via Shapiro–Wilk and Jarque–Bera tests. Quantile-quantile (Q-Q) plots visually confirmed approximate normality, supporting parametric inference.

---

## SI-1.4 Notation, Symbols, and Conventions

**Primary symbols:**

*Λ³ thermal-softening model (SI-2):*
- E(T): Young's modulus at temperature T [Pa]
- E₀: Young's modulus at reference temperature T₀ = 300 K [Pa]
- λ_base: Baseline harmonic attenuation coefficient [dimensionless]
- κ: Anharmonicity parameter (softening acceleration) [dimensionless]
- α: Linear thermal expansion coefficient [K⁻¹]
- T_m: Melting temperature [K]
- ΔT: Temperature difference T − T₀ [K]

*Elastic moduli (SI-2, SI-4):*
- G(T): Shear modulus at temperature T [Pa]
- K(T): Bulk modulus at temperature T [Pa]
- ν: Poisson's ratio [dimensionless]
- f_G: Shear retention factor G(T)/G(300 K) [dimensionless]

*Lindemann melting criterion (SI-3, SI-4):*
- δ_L: Lindemann ratio ⟨u²⟩^(1/2) / r_nn [dimensionless]
- ⟨u²⟩: Mean-square atomic displacement [m²]
- r_nn: Nearest-neighbor interatomic distance [m]
- T_m: Melting point [K]

*Debye–Waller theory (SI-4):*
- k_B: Boltzmann constant = 1.38065×10⁻²³ J/K
- M: Atomic mass [kg]
- ω: Phonon angular frequency [rad/s]
- v_t: Transverse (shear) sound velocity [m/s]
- v_l: Longitudinal sound velocity [m/s]
- k_D: Debye wavevector [m⁻¹]
- n: Atomic number density [m⁻³]

*Energy-density ratio (SI-5):*
- EDR: Energy-density ratio k_B T / (K V_atom) [dimensionless]
- V_atom: Atomic volume [m³]
- Λ: Dimensionless energy-density ratio K/|V| [dimensionless]

*Gibbs–Thomson nanoparticle melting (SI-6):*
- r: Nanoparticle radius [m or nm]
- T_m(r): Size-dependent melting point [K]
- T_m^∞: Bulk melting point [K]
- ΔT_m: Melting-point depression T_m^∞ − T_m(r) [K]
- a, b: Capillary length parameters (1/r and 1/r² terms) [m]
- Z: Coordination number [dimensionless]

**Crystal structure notation:**
- BCC: Body-centered cubic (e.g., Fe, W)
- FCC: Face-centered cubic (e.g., Cu, Al, Ni)
- HCP: Hexagonal close-packed (e.g., Ti, Mg)

**Statistical notation:**
- R²: Coefficient of determination (uncentered for through-origin models) [dimensionless]
- RMSE: Root-mean-square error [units of dependent variable]
- CI: Confidence interval (95% unless stated otherwise)
- SE: Standard error
- LOO: Leave-one-out cross-validation

**Unit conventions:**
All physical quantities are reported in SI units unless explicitly stated. Temperature is given in Kelvin [K], with Celsius [°C] provided for reference where appropriate. Atomic-scale lengths (lattice constants, interatomic distances) are reported in Ångströms [Å] for convenience, with conversions to meters [m] used in calculations.

---

## SI-1.5 Data Quality Assessment Criteria

Experimental and computational data sources were classified into two tiers based on measurement reliability, reproducibility, and consistency across independent studies:

**Tier 1 (High confidence):**
- Data from peer-reviewed experimental studies with explicit error bars
- Multiple independent measurements showing consistency (< 10% variation)
- Well-documented measurement techniques (ultrasonic pulse-echo, resonant ultrasound spectroscopy, X-ray/neutron diffraction, differential scanning calorimetry)
- High-purity materials (≥ 99.9% purity for elemental metals)
- Temperature-controlled environments with minimal oxidation/contamination

*Examples:* High-temperature elastic moduli from ultrasonic methods (Ledbetter series, Ogi et al.); nanoparticle melting from in-situ TEM with heating stage (Buffat & Borel, Zhang et al.); ab initio DFT calculations for Lindemann parameters (Kavarnadze et al., Vočadlo & Alfè).

**Tier 2 (Moderate confidence):**
- Single-source data without independent validation
- Manufacturer datasheets or engineering handbooks (ASM, NIST) where original measurement conditions are not fully specified
- Extrapolated or interpolated values from limited temperature ranges
- Molecular dynamics (MD) simulations using empirical potentials (EAM, Finnis–Sinclair) without DFT validation

**Data exclusion criteria:**
Data points were excluded if:
- Measurements were performed in phase-transition regions (e.g., α→γ transition in Fe at 1185 K, magnetic Curie points)
- Evident sample degradation (oxidation, grain growth, phase separation) during high-temperature measurements
- Inconsistency exceeding ±20% compared to multiple independent sources
- Nanoparticle data below r < 2 nm, where quantum confinement and non-bulk behavior dominate

**Literature selection for Lindemann parameters (SI-3):**
Priority was given to:
1. Recent ab initio calculations (post-2000) using DFT-based phonon theory
2. Experimental determinations via neutron/X-ray scattering (Debye–Waller factors)
3. Systematic theoretical studies with explicit structure-dependent analysis (Shapiro 1970, Kuramoto et al. 2010)

Original literature reporting δ_L relative to lattice constant a was converted to nearest-neighbor distance r_nn for consistency (e.g., FCC: r_nn = a/√2, conversion factor √2 ≈ 1.414).

---

## SI-1.6 Code and Data Availability Statement

**Source code:**
All computational scripts for data processing, Λ³ fitting, Lindemann predictions, EDR regression, and figure generation are available at:

GitHub repository: `https://github.com/[username]/lambda3-melting-criterion`

The repository includes:
- Python scripts for all calculations
  - `lindemann_predictor.py` (SI-4: Born-struct predictions)
  - `gibbs_thomson_analysis.py` (SI-6: Nanoparticle analysis)
  - `lambda3_thermal_softening.py` (SI-2: Thermal softening fits)
- Requirements file (`requirements.txt`) specifying dependencies
- README with installation and usage instructions

**Data files:**
All experimental and computational data used in this study are provided as CSV files in the repository `data/` directory:

1. `thermal_properties.csv` – Temperature-dependent E(T), σ_y(T), α(T) for 7 metals (SI-2)
2. `material_summary.csv` – Λ³ fitted parameters (λ_base, κ) and material properties (SI-2)
3. `lindemann_parameters_literature.csv` – Experimental δ_L values with literature sources (SI-3)
4. `lindemann_physical_results.csv` – Born-struct predictions and validation (SI-4)
5. `nanoparticle_melting_gibbs_thomson.csv` – Size-dependent melting data for 7 materials (SI-6)

**Persistent archival:**
Upon acceptance, a snapshot of the repository will be archived on Zenodo with a citable DOI.

**Correspondence:**
Requests for additional data or clarification should be directed to:
- Masamichi Iizumi (Miosync Inc.): [email to be added]

**License:**
All code and data are released under the MIT License, permitting unrestricted use with attribution.

---

## End of SI-1

Subsequent sections (SI-2 through SI-6) provide detailed technical documentation as outlined in Section SI-1.1.
