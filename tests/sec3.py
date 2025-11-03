import numpy as np
from scipy import interpolate
from scipy.stats import poisson
import matplotlib.pyplot as plt

from eemod import maximum_proton_energy, neutrino_spectra

GeV_to_erg = 1.602e-3  # 1 GeV = 1.602e-3 erg
eV_to_erg = 1.60217663e-12

Gamma = 30.0
r_diss = 1.0e14  # cm
B = 2.920420e3   # G
alpha = 0.5
beta = 2.0
eps_pk_comoving_erg = 5.341e-11  
eps_pk_eV = eps_pk_comoving_erg / eV_to_erg
A = 1.73920e+14 / eV_to_erg  

E_gamma_iso = 3.800e51  # erg
xi_p = 10.0
d_L_Mpc = 300.0


proton_results = maximum_proton_energy(
    Gamma=Gamma,
    r_diss=r_diss,
    B=B,
    eps_pk_eV=eps_pk_eV,
    alpha=alpha,
    beta=beta,
    A=A,
    eta_acc=1.0,
    n_points=5000
)

# Calculate neutrino spectra
neu_results = neutrino_spectra(
    proton_results,
    E_gamma_iso=E_gamma_iso,
    xi_p=xi_p,
    d_L_Mpc=d_L_Mpc
)


E_nu_obs_erg = neu_results['E_nu_obs']  # [erg]
phi_numu_earth = neu_results['phi_numu_earth']  # [#/(erg·cm²)]

E_nu_GeV = E_nu_obs_erg / GeV_to_erg

phi_numu_per_GeV = phi_numu_earth * GeV_to_erg


aeff_filename = "Aeff_all.txt"

try:
    data_aeff = np.loadtxt(aeff_filename)
    E_aeff_GeV = data_aeff[:, 0]
    Aeff_cm2   = data_aeff[:, 1]
    print(f"✓ Loaded {aeff_filename}")
except Exception as e:
    print(f"❌ ERROR loading {aeff_filename}: {e}")
    exit(1)




E_min = max(E_nu_GeV.min(), E_aeff_GeV.min(), 100)
E_max = min(E_nu_GeV.max(), E_aeff_GeV.max())


E_common_GeV = np.logspace(np.log10(E_min), np.log10(E_max), 500)

interp_phi = interpolate.interp1d(
    np.log10(E_nu_GeV),
    np.log10(phi_numu_per_GeV + 1e-100),
    kind='linear',
    bounds_error=False,
    fill_value=-np.inf
)
phi_common = 10**interp_phi(np.log10(E_common_GeV))
phi_common = np.maximum(phi_common, 0)

interp_aeff = interpolate.interp1d(
    np.log10(E_aeff_GeV),
    np.log10(Aeff_cm2 + 1e-10),
    kind='linear',
    bounds_error=False,
    fill_value=-np.inf
)
Aeff_common = 10**interp_aeff(np.log10(E_common_GeV))
Aeff_common = np.maximum(Aeff_common, 0)


integrand = phi_common * Aeff_common

N_mu = np.trapz(integrand, E_common_GeV)


print(f"\n{'='*80}")
print(f"RESULTS: EXPECTED NUMBER OF EVENTS")
print(f"{'='*80}")
print(f"\nAt d_L = {d_L_Mpc} Mpc:")
print(f"  N_μ = {N_mu:.6f}")



P_geq_1 = 1 - poisson.cdf(0, mu=N_mu)
P_geq_2 = 1 - poisson.cdf(1, mu=N_mu)

print(f"\n  P(N_μ ≥ 1) = {P_geq_1:.6f} ({P_geq_1*100:.3f}%)")
print(f"  P(N_μ ≥ 2) = {P_geq_2:.6f} ({P_geq_2*100:.3f}%)")

print(f"\nPoisson distribution:")
print(f"  {'k':<3} {'P(N_μ = k)':<15} {'P(N_μ ≥ k)':<15}")
print(f"  {'-'*35}")
for k in range(6):
    p_k = poisson.pmf(k, mu=N_mu)
    p_geq_k = 1 - poisson.cdf(k-1, mu=N_mu)
    print(f"  {k:<3} {p_k:<15.6f} {p_geq_k:<15.6f}")



d_L_array = np.array([50, 100, 150, 200, 300, 400, 500, 600, 800, 1000])
N_mu_array = N_mu * (d_L_Mpc / d_L_array)**2
P_1_array = 1 - poisson.cdf(0, mu=N_mu_array)
P_2_array = 1 - poisson.cdf(1, mu=N_mu_array)

print(f"\n  d_L [Mpc]  |  N_μ        |  P(N_μ≥1)   |  P(N_μ≥2)")
print(f"  {'-'*55}")
for i, d_L in enumerate(d_L_array):
    print(f"  {d_L:4.0f}       |  {N_mu_array[i]:8.5f}   |  {P_1_array[i]:8.5f}  |  {P_2_array[i]:8.5f}")

# =============================================================================
# STEP 7: Comparison with Paper
# =============================================================================
print(f"\n{'='*80}")
print(f"COMPARISON WITH PAPER (Table 2)")
print(f"{'='*80}")

print(f"""
EE-mod (single Γ = 30)
  N_μ = {N_mu:.6f}
  P(N_μ ≥ 1) = {P_geq_1:.4f} ({P_geq_1*100:.2f}%)

Paper expectations (Table 2, d_L = 300 Mpc):
  EE-mod (single Γ=30): P(N_μ≥1) ≈ 0.04 (4%)
  EE-mod-dist-A (σ_Γ=2): P(N_μ≥1) ≈ 0.07 (7%)
  EE-mod-dist-B (σ_Γ=4): P(N_μ≥1) ≈ 0.11 (11%)
""")
