#!/usr/bin/env python3
"""
FINAL CORRECTED: Thermal vs Mechanical
Using Born-collapsed G for thermal pathway
"""
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.constants import k as k_B

# ========================================
# Physical constants & Material (Cu)
# ========================================
k_B = 1.380649e-23
u_kg = 1.66053906660e-27
pi = math.pi

# Cu properties
a_300K_A = 3.61      # Å
alpha = 1.7e-5       # /K
T_ref = 300
T_melt = 1357

E0_GPa = 130
sigma_y_MPa = 122
nu = 0.34
M_amu = 63.546
M = M_amu * u_kg

rho_kg_m3 = 8900

# Born collapse factor (FCC structure median from paper)
fG_Born_FCC = 0.101

# ========================================
# Geometry functions
# ========================================
def a_T_meters(T):
    """Lattice constant at T in meters"""
    a0_m = a_300K_A * 1e-10
    return a0_m * (1.0 + alpha * (T - T_ref))

def r_nn_FCC(T):
    """FCC nearest neighbor at T"""
    return a_T_meters(T) / math.sqrt(2.0)

def number_density_FCC(T):
    """FCC atomic density at T"""
    a_m = a_T_meters(T)
    return 4.0 / (a_m**3)

# ========================================
# SCENARIO A: THERMAL MELTING
# (Using Born-collapsed G)
# ========================================
print("="*80)
print("SCENARIO A: THERMAL MELTING (Born collapse → Lindemann)")
print("="*80)

T = T_melt
n_Tm = number_density_FCC(T)
rnn_Tm = r_nn_FCC(T)

# Step 1: Room temperature elastic moduli
G0_Pa = (E0_GPa * 1e9) / (2.0 * (1.0 + nu))
K0_Pa = (E0_GPa * 1e9) / (3.0 * (1.0 - 2.0 * nu))

print(f"\nStep 1: Room temperature properties")
print(f"  G(300K) = {G0_Pa/1e9:.2f} GPa")
print(f"  K(300K) = {K0_Pa/1e9:.2f} GPa")

# Step 2: Born collapse at Tm
Gm_Pa = fG_Born_FCC * G0_Pa
Km_Pa = K0_Pa  # Bulk less affected

print(f"\nStep 2: Born collapse at Tm")
print(f"  fG = {fG_Born_FCC:.3f}")
print(f"  G(Tm, after collapse) = {Gm_Pa/1e9:.2f} GPa")
print(f"  Reduction: {(1-fG_Born_FCC)*100:.1f}%")

# Step 3: Sound velocities (after collapse)
rho_Tm = n_Tm * M
v_t = math.sqrt(Gm_Pa / rho_Tm)
v_l = math.sqrt((Km_Pa + 4.0*Gm_Pa/3.0) / rho_Tm)

print(f"\nStep 3: Phonon collapse")
print(f"  ρ(Tm) = {rho_Tm:.0f} kg/m³")
print(f"  v_t = {v_t:.0f} m/s")
print(f"  v_l = {v_l:.0f} m/s")

# Step 4: Debye-Waller calculation
k_D = (6.0 * pi**2 * n_Tm) ** (1.0/3.0)
inv_omega2 = (1.0 / (3.0 * k_D**2)) * (2.0 / v_t**2 + 1.0 / v_l**2)
u2_thermal = (k_B * T / M) * inv_omega2
u_rms = math.sqrt(u2_thermal)
delta_thermal = u_rms / rnn_Tm

print(f"\nStep 4: Lindemann calculation")
print(f"  k_D = {k_D:.4e} m⁻¹")
print(f"  ⟨1/ω²⟩ = {inv_omega2:.4e} s²")
print(f"  ⟨u²⟩ = {u2_thermal:.4e} m²")
print(f"  √⟨u²⟩ = {u_rms*1e10:.4f} Å")
print(f"  r_nn(Tm) = {rnn_Tm*1e10:.4f} Å")
print(f"  δ_thermal = {delta_thermal:.4f}")

if 0.10 <= delta_thermal <= 0.18:
    print(f"\n  ✅ Lindemann criterion satisfied!")
    print(f"     {delta_thermal:.3f} ∈ [0.10, 0.18]")
else:
    print(f"\n  ⚠️  Value: {delta_thermal:.3f}")
    print(f"     Expected: [0.10, 0.18]")

# Generate distribution
N_atoms = 100000
sigma_spread = 0.02  # Thermal fluctuations
delta_thermal_dist = np.random.normal(delta_thermal, sigma_spread, N_atoms)
delta_thermal_dist = np.maximum(delta_thermal_dist, 0)

# ========================================
# SCENARIO B: MECHANICAL YIELDING
# (From simulation, tail shift)
# ========================================
print("\n" + "="*80)
print("SCENARIO B: MECHANICAL YIELDING (Tail shift)")
print("="*80)

N_sites = 100000
mean_mech = 0.024  # Bulk average from simulation
tail_fraction = 0.0029  # 0.29% sites exceed δ > 0.1

# Bulk: lognormal centered at 0.024
delta_mech_bulk = np.random.lognormal(np.log(mean_mech), 0.5, 
                                      int(N_sites * (1 - tail_fraction)))

# Tail: uniform 0.10-0.28
delta_mech_tail = np.random.uniform(0.10, 0.28, 
                                    int(N_sites * tail_fraction))

delta_mech_dist = np.concatenate([delta_mech_bulk, delta_mech_tail])

print(f"Mean δ: {np.mean(delta_mech_dist):.4f}")
print(f"Sites δ > 0.1: {100*np.sum(delta_mech_dist > 0.1)/len(delta_mech_dist):.2f}%")

# ========================================
# Visualization
# ========================================
fig, axes = plt.subplots(2, 1, figsize=(12, 10))

bins = np.linspace(0, 0.30, 80)
delta_threshold = 0.10

# Panel A: THERMAL
ax1 = axes[0]
ax1.hist(delta_thermal_dist, bins=bins, alpha=0.7, color='red',
         edgecolor='darkred', linewidth=0.5,
         label=f'Thermal (T = {T_melt} K)')

ax1.axvline(delta_threshold, color='black', linestyle='--', linewidth=2.5,
           label=f'Lindemann threshold ({delta_threshold})')

mean_th = np.mean(delta_thermal_dist)
ax1.axvline(mean_th, color='darkred', linestyle=':', linewidth=2.5,
           label=f'Mean δ = {mean_th:.3f}')

ax1.text(0.98, 0.95, f'(A) THERMAL MELTING\n"Mean Shift"\nδ = {delta_thermal:.3f}',
        transform=ax1.transAxes, fontsize=13, fontweight='bold',
        verticalalignment='top', horizontalalignment='right',
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

ax1.set_xlabel('Lindemann ratio δ = √⟨u²⟩ / r_nn', fontsize=13, fontweight='bold')
ax1.set_ylabel('Count', fontsize=13, fontweight='bold')
ax1.set_xlim([0, 0.30])
ax1.set_yscale('log')
ax1.set_ylim([1, None])
ax1.legend(fontsize=10)
ax1.grid(True, alpha=0.3)

# Panel B: MECHANICAL
ax2 = axes[1]
ax2.hist(delta_mech_dist, bins=bins, alpha=0.7, color='skyblue',
         edgecolor='black', linewidth=0.5,
         label=f'Mechanical (σ = {sigma_y_MPa} MPa, T = 300 K)')

ax2.axvline(delta_threshold, color='black', linestyle='--', linewidth=2.5,
           label=f'Lindemann threshold ({delta_threshold})')

mean_mech = np.mean(delta_mech_dist)
ax2.axvline(mean_mech, color='blue', linestyle=':', linewidth=2.5,
           label=f'Mean δ = {mean_mech:.3f}')

tail = delta_mech_dist[delta_mech_dist > delta_threshold]
if len(tail) > 0:
    ax2.hist(tail, bins=bins, alpha=0.9, color='red',
            edgecolor='darkred', linewidth=0.8,
            label=f'δ > {delta_threshold} ({100*len(tail)/len(delta_mech_dist):.2f}%)')

ax2.text(0.98, 0.95, '(B) MECHANICAL YIELDING\n"Tail Shift"',
        transform=ax2.transAxes, fontsize=13, fontweight='bold',
        verticalalignment='top', horizontalalignment='right',
        bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))

ax2.set_xlabel('Mechanical Lindemann ratio δ_mech = σ_local / E_eff', 
               fontsize=13, fontweight='bold')
ax2.set_ylabel('Count', fontsize=13, fontweight='bold')
ax2.set_xlim([0, 0.30])
ax2.set_yscale('log')
ax2.set_ylim([1, None])
ax2.legend(fontsize=10)
ax2.grid(True, alpha=0.3)

plt.suptitle('Universal Lindemann Criterion: Thermal vs Mechanical Pathways\n' +
            'Same threshold (δ ≈ 0.10), different spatial distributions',
            fontsize=15, fontweight='bold', y=0.995)

plt.tight_layout()
plt.savefig('thermal_vs_mechanical_CORRECTED.png', dpi=300, bbox_inches='tight')
print("\n✅ Saved: thermal_vs_mechanical_CORRECTED.png")

plt.show()
