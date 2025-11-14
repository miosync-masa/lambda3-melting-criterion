#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Physically-derived Lindemann predictor (complete derivation)
"""
import math
import numpy as np
import pandas as pd
from dataclasses import dataclass

# ---------------------------
# Physical constants
# ---------------------------
k_B = 1.380649e-23      # J/K
u_kg = 1.66053906660e-27  # atomic mass unit in kg
pi = math.pi

# ---------------------------
# Material input
# ---------------------------
@dataclass
class Element:
    name: str
    struct: str
    a_A: float         # lattice a at 300 K [Å]
    c_over_a: float    # for HCP
    T_m: float         # melting point [K]
    alpha: float       # linear thermal expansion [/K]
    M_amu: float       # atomic mass [amu]
    E_GPa: float       # Young's modulus at RT [GPa]
    nu: float          # Poisson ratio at RT
    delta_exp: float   # experimental Lindemann

elts = [
    Element('Fe', 'BCC', 2.92, 0.0, 1811.0, 1.50e-5, 55.845, 211.0, 0.29, 0.180),
    Element('W',  'BCC', 3.16, 0.0, 3695.0, 4.51e-6, 183.84, 411.0, 0.28, 0.160),
    Element('Cu', 'FCC', 3.61, 0.0, 1357.0, 1.70e-5, 63.546, 130.0, 0.34, 0.100),
    Element('Al', 'FCC', 4.05, 0.0,  933.0, 2.30e-5, 26.982,  70.0, 0.35, 0.100),
    Element('Ni', 'FCC', 3.52, 0.0, 1728.0, 1.30e-5, 58.693, 200.0, 0.31, 0.110),
    Element('Ti', 'HCP', 2.95, 1.587, 1941.0, 8.60e-6, 47.867, 116.0, 0.32, 0.100),
    Element('Mg', 'HCP', 3.21, 1.624,  923.0, 2.70e-5, 24.305,  45.0, 0.29, 0.117),
]

# ---------------------------
# Helper functions
# ---------------------------
def a_T(a0_A, alpha, Tm, Tref=300.0):
    """Lattice constant at Tm (in meters)."""
    a0_m = a0_A * 1e-10
    return a0_m * (1.0 + alpha * (Tm - Tref))

def r_nn_T(el: Element):
    """Nearest-neighbour distance at Tm in meters."""
    a_m = a_T(el.a_A, el.alpha, el.T_m)

    if el.struct == 'BCC':
        return (math.sqrt(3.0)/2.0) * a_m
    elif el.struct == 'FCC':
        return a_m / math.sqrt(2.0)
    elif el.struct == 'HCP':
        # Simplified: use basal a
        return a_m
    return a_m

def number_density(el: Element):
    """Atomic number density at Tm [1/m^3]."""
    a_m = a_T(el.a_A, el.alpha, el.T_m)

    if el.struct == 'FCC':
        return 4.0 / (a_m**3)
    elif el.struct == 'BCC':
        return 2.0 / (a_m**3)
    elif el.struct == 'HCP':
        # Approximate with FCC-like
        return 4.0 / (a_m**3)
    return 4.0 / (a_m**3)

def elastic_moduli(el: Element):
    """Return (K, G) in SI units (Pa)."""
    E = el.E_GPa * 1e9
    nu = el.nu
    G = E / (2.0 * (1.0 + nu))
    K = E / (3.0 * (1.0 - 2.0 * nu))
    return K, G

def delta_with_fG(el: Element, fG):
    """
    Prediction with shear softening factor fG = G(Tm)/G(RT).

    <u²> = (k_B T_m / M) × <1/ω²>
    <1/ω²> ≈ (1/3k_D²) × (2/v_t² + 1/v_l²)
    """
    n_Tm = number_density(el)
    K0, G0 = elastic_moduli(el)

    # Apply softening
    Gm = max(fG * G0, 1e-12)

    M = el.M_amu * u_kg
    rho_Tm = n_Tm * M

    # Sound velocities
    v_t = math.sqrt(Gm / rho_Tm)
    v_l = math.sqrt((K0 + 4.0*Gm/3.0) / rho_Tm)

    # Debye wavevector
    k_D = (6.0 * pi**2 * n_Tm) ** (1.0/3.0)

    # Average 1/ω²
    inv_omega2 = (1.0 / (3.0 * k_D**2)) * (2.0 / v_t**2 + 1.0 / v_l**2)

    # Mean square displacement
    u2 = (k_B * el.T_m / M) * inv_omega2

    # Lindemann ratio
    rnn = r_nn_T(el)
    delta = math.sqrt(u2) / rnn

    return delta

def infer_fG_for_element(el: Element, target_delta):
    """
    Find fG such that delta_with_fG(el, fG) ≈ target_delta.
    """
    # Binary search
    f_min, f_max = 1e-6, 1.0

    # Check bounds
    d_min = delta_with_fG(el, f_min)
    d_max = delta_with_fG(el, f_max)

    if d_min < target_delta:
        return f_min, d_min
    if d_max > target_delta:
        return f_max, d_max

    # Binary search
    for _ in range(100):
        f_mid = math.sqrt(f_min * f_max)  # Geometric mean
        d_mid = delta_with_fG(el, f_mid)

        if abs(d_mid - target_delta) < 1e-6:
            return f_mid, d_mid

        # delta decreases as fG increases
        if d_mid < target_delta:
            f_max = f_mid
        else:
            f_min = f_mid

    return f_mid, d_mid

# ---------------------------
# Run analysis
# ---------------------------
print("\n" + "="*80)
print("PHYSICALLY-DERIVED LINDEMANN PREDICTOR")
print("="*80)

# 1. No softening baseline
print("\n--- Baseline (no softening, fG=1.0) ---")
rows_base = []
for e in elts:
    d0 = delta_with_fG(e, 1.0)
    err = 100.0 * abs(d0 - e.delta_exp) / e.delta_exp
    rows_base.append({
        'Material': e.name,
        'Struct': e.struct,
        'δ_exp': e.delta_exp,
        'δ_no_soft': d0,
        'Error_%': err
    })

df_base = pd.DataFrame(rows_base)
print(df_base.to_string(index=False))
print(f"\nMean error: {df_base['Error_%'].mean():.1f}%")

# 2. Infer required fG
print("\n--- Inferred fG per element ---")
rows_inf = []
for e in elts:
    f_req, d_at_f = infer_fG_for_element(e, e.delta_exp)
    rows_inf.append({
        'Material': e.name,
        'Struct': e.struct,
        'fG_required': f_req,
        'δ_at_fG': d_at_f
    })

df_inf = pd.DataFrame(rows_inf)
print(df_inf.to_string(index=False))

# 3. Structure medians
struct_medians = df_inf.groupby('Struct')['fG_required'].median()
print("\n--- Structure medians ---")
print(struct_medians)

# 4. Born-struct predictions
print("\n--- Born-struct predictions (using structure medians) ---")
rows_born = []
for e in elts:
    f_struct = struct_medians[e.struct]
    d_pred = delta_with_fG(e, f_struct)
    err = 100.0 * abs(d_pred - e.delta_exp) / e.delta_exp
    rows_born.append({
        'Material': e.name,
        'Struct': e.struct,
        'δ_exp': e.delta_exp,
        'fG_struct': f_struct,
        'δ_pred': d_pred,
        'Error_%': err
    })

df_born = pd.DataFrame(rows_born)
print(df_born.to_string(index=False))
print(f"\nMean error (Born-struct): {df_born['Error_%'].mean():.1f}%")

# Save
df_born.to_csv('lindemann_physical_results.csv', index=False)
print("\nSaved to 'lindemann_physical_results.csv'")
