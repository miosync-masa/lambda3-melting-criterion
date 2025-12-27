#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Z³ scaling law analysis for shear retention factor f_G
Hypothesis: f_G ∝ (Z/12)³

Testing whether coordination number explains structure-dependent collapse
"""
import numpy as np
import pandas as pd

# =============================================================================
# Experimental data from paper
# =============================================================================
data = {
    'Material': ['Fe', 'W', 'Cu', 'Al', 'Ni', 'Ti', 'Mg'],
    'Struct':   ['BCC', 'BCC', 'FCC', 'FCC', 'FCC', 'HCP', 'HCP'],
    'Z':        [8, 8, 12, 12, 12, 12, 12],  # Nominal coordination
    'c_over_a': [None, None, None, None, None, 1.587, 1.624],
    'f_G_req':  [0.026, 0.028, 0.111, 0.101, 0.079, 0.081, 0.077],  # From paper
    'delta_exp': [0.180, 0.160, 0.100, 0.100, 0.110, 0.100, 0.117],
}

df = pd.DataFrame(data)

print("="*70)
print("Z³ SCALING LAW ANALYSIS")
print("="*70)

# =============================================================================
# Step 1: Check Z³ scaling for BCC vs FCC
# =============================================================================
print("\n--- Step 1: BCC vs FCC comparison ---")

f_G_BCC_mean = df[df['Struct'] == 'BCC']['f_G_req'].mean()
f_G_FCC_mean = df[df['Struct'] == 'FCC']['f_G_req'].mean()
f_G_HCP_mean = df[df['Struct'] == 'HCP']['f_G_req'].mean()

print(f"f_G mean (BCC): {f_G_BCC_mean:.4f}")
print(f"f_G mean (FCC): {f_G_FCC_mean:.4f}")
print(f"f_G mean (HCP): {f_G_HCP_mean:.4f}")

# Z³ prediction
Z_BCC, Z_FCC = 8, 12
ratio_Z3 = (Z_BCC / Z_FCC) ** 3
ratio_actual = f_G_BCC_mean / f_G_FCC_mean

print(f"\n(Z_BCC/Z_FCC)³ = ({Z_BCC}/{Z_FCC})³ = {ratio_Z3:.4f}")
print(f"f_G_BCC / f_G_FCC (actual) = {ratio_actual:.4f}")
print(f"Error: {100*abs(ratio_Z3 - ratio_actual)/ratio_actual:.1f}%")

# =============================================================================
# Step 2: If Z³ holds, what is Z_eff for HCP?
# =============================================================================
print("\n--- Step 2: Infer Z_eff for HCP ---")

# Using FCC as reference: f_G ∝ (Z/12)³
# f_G_HCP / f_G_FCC = (Z_eff / 12)³
# Z_eff = 12 × (f_G_HCP / f_G_FCC)^(1/3)

ratio_HCP_FCC = f_G_HCP_mean / f_G_FCC_mean
Z_eff_HCP = 12 * (ratio_HCP_FCC) ** (1/3)

print(f"f_G_HCP / f_G_FCC = {ratio_HCP_FCC:.4f}")
print(f"Z_eff (HCP) = 12 × {ratio_HCP_FCC:.4f}^(1/3) = {Z_eff_HCP:.2f}")
print(f"Nominal Z (HCP) = 12")
print(f"Effective reduction: {12 - Z_eff_HCP:.2f} ({100*(12-Z_eff_HCP)/12:.1f}%)")

# =============================================================================
# Step 3: Physical interpretation - c/a ratio effect?
# =============================================================================
print("\n--- Step 3: c/a ratio analysis for HCP ---")

# Ideal HCP: c/a = sqrt(8/3) ≈ 1.633
c_a_ideal = np.sqrt(8/3)
print(f"Ideal c/a ratio: {c_a_ideal:.3f}")

hcp_data = df[df['Struct'] == 'HCP'].copy()
for _, row in hcp_data.iterrows():
    c_a = row['c_over_a']
    deviation = (c_a - c_a_ideal) / c_a_ideal * 100
    print(f"{row['Material']}: c/a = {c_a:.3f}, deviation from ideal = {deviation:+.1f}%")

# Hypothesis: c/a deviation weakens c-axis bonds
# Z_eff = 6 (basal) + 6 × g(c/a) where g < 1 for non-ideal
print("\n--- Hypothesis: Z_eff = 6 + 6 × g(c/a) ---")

for _, row in hcp_data.iterrows():
    c_a = row['c_over_a']
    f_G = row['f_G_req']
    
    # Infer Z_eff from f_G
    Z_eff_ind = 12 * (f_G / f_G_FCC_mean) ** (1/3)
    
    # If Z_eff = 6 + 6*g, then g = (Z_eff - 6)/6
    g = (Z_eff_ind - 6) / 6
    
    print(f"{row['Material']}: c/a={c_a:.3f}, f_G={f_G:.3f}, Z_eff={Z_eff_ind:.2f}, g={g:.3f}")

# =============================================================================
# Step 4: Universal Z³ model with Z_eff
# =============================================================================
print("\n--- Step 4: Universal Z³ model predictions ---")

# Define Z_eff values
Z_eff_dict = {
    'BCC': 8,
    'FCC': 12,
    'HCP': Z_eff_HCP  # Use average inferred value
}

# Reference: FCC with Z=12
f_G_ref = f_G_FCC_mean
Z_ref = 12

print(f"\nReference: f_G(FCC, Z=12) = {f_G_ref:.4f}")
print(f"\nZ³ scaling law: f_G = {f_G_ref:.4f} × (Z_eff/12)³")

print("\nPredicted f_G values:")
for struct in ['BCC', 'FCC', 'HCP']:
    Z_eff = Z_eff_dict[struct]
    f_G_pred = f_G_ref * (Z_eff / Z_ref) ** 3
    f_G_actual = df[df['Struct'] == struct]['f_G_req'].mean()
    err = 100 * abs(f_G_pred - f_G_actual) / f_G_actual
    print(f"  {struct}: Z_eff={Z_eff:.1f}, f_G_pred={f_G_pred:.4f}, f_G_actual={f_G_actual:.4f}, error={err:.1f}%")

# =============================================================================
# Step 5: Per-element validation
# =============================================================================
print("\n--- Step 5: Per-element validation ---")

results = []
for _, row in df.iterrows():
    Z_eff = Z_eff_dict[row['Struct']]
    f_G_pred = f_G_ref * (Z_eff / Z_ref) ** 3
    f_G_actual = row['f_G_req']
    err = 100 * abs(f_G_pred - f_G_actual) / f_G_actual
    
    results.append({
        'Material': row['Material'],
        'Struct': row['Struct'],
        'Z_eff': Z_eff,
        'f_G_pred': f_G_pred,
        'f_G_actual': f_G_actual,
        'Error_%': err
    })

df_results = pd.DataFrame(results)
print(df_results.to_string(index=False))
print(f"\nMean absolute error: {df_results['Error_%'].mean():.1f}%")

# =============================================================================
# Step 6: Physical formula summary
# =============================================================================
print("\n" + "="*70)
print("SUMMARY: PHYSICAL DERIVATION OF f_G")
print("="*70)

print("""
Proposed scaling law:

    f_G = f_G⁰ × (Z_eff / 12)³

where:
    f_G⁰ = 0.097 (FCC reference, close-packed limit)
    
    Z_eff = {
        8    for BCC (open structure)
        12   for FCC (close-packed cubic)
        ~11  for HCP (close-packed hexagonal with c/a anisotropy)
    }

Physical interpretation:
    
    1st power (Z¹): Bond energy density ∝ Z
    2nd power (Z²): Angular constraints (shear resistance) ∝ Z²  
    3rd power (Z³): Cooperative collapse probability (rigidity percolation)
    
    Combined: f_G ∝ Z³ (3D network rigidity)
""")

print("="*70)
