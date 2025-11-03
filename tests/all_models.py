import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from eemod import maximum_proton_energy, neutrino_spectra


c = 2.99792458e10
m_e = 9.1093897e-28
m_pion = 2.48832e-25
m_muon = 1.88353e-25
sigma_T = 6.6524587e-25
tau_pion = 2.6e-8
tau_muon = 2.197e-6
eV_to_erg = 1.602e-12


models = {
    "EE-mod":   {"Gamma": 30,  "L_obs": 3e48, "E_obs": 1e51, "r_diss": 1e14,  "Epk_keV": 1.0},
    "EE-opt":   {"Gamma": 10,  "L_obs": 3e48, "E_obs": 1e51, "r_diss": 3e13,  "Epk_keV": 10.0},
    "Prompt":   {"Gamma": 1e3, "L_obs": 1e51, "E_obs": 1e51, "r_diss": 3e13,  "Epk_keV": 500.0},
    "Flare":    {"Gamma": 30,  "L_obs": 1e48, "E_obs": 3e50, "r_diss": 3e14,  "Epk_keV": 0.3},
    "Plateau":  {"Gamma": 30,  "L_obs": 1e47, "E_obs": 3e50, "r_diss": 3e14,  "Epk_keV": 0.1},
}

def run_model(name, params, alpha=0.5, beta=2.0, xi_B=0.1, xi_p=10.0, d_L_Mpc=300.0, eta_acc=1.0):
    Gamma, r_diss, L_obs, E_gamma_iso_star, Epk_keV = (
        params["Gamma"], params["r_diss"], params["L_obs"], params["E_obs"], params["Epk_keV"]
    )
    eps_pk_eV = (Epk_keV * 1e3) / Gamma


    def dn_depsilon(eps, A, eps_pk, alpha, beta):
        return A * np.where(eps < eps_pk, (eps / eps_pk) ** (-alpha), (eps / eps_pk) ** (-beta))

    def U(A, eps_pk, alpha, beta, eps_min, eps_max):
        integrand = lambda e: e * dn_depsilon(e, A, eps_pk, alpha, beta)
        res, _ = integrate.quad(integrand, eps_min, eps_max, limit=200)
        return res
  
    if name == "Prompt":
       eps_min_obs = (10 / Gamma) * 1000
       eps_max_obs = (1e3 / Gamma) * 1000
    else:
        eps_min_obs =  (0.3 / Gamma) * 1000
        eps_max_obs =  (10.0 / Gamma) * 1000

    eps_min_full = 0.1
    eps_max_full = 1e6
    U_obs_norm = U(1.0, eps_pk_eV, alpha, beta, eps_min_obs, eps_max_obs)
    U_full_norm = U(1.0, eps_pk_eV, alpha, beta, eps_min_full, eps_max_full)
    corr = U_full_norm / U_obs_norm

    coeff = 4 * np.pi * c * Gamma**2 * r_diss**2 * U_obs_norm * eV_to_erg
    A = L_obs / coeff / eV_to_erg
    U_full_actual = A * U_full_norm * eV_to_erg
    L_iso = 4 * np.pi * c * Gamma**2 * r_diss**2 * U_full_actual * eV_to_erg
    B = np.sqrt(2 * L_iso * xi_B / (c * Gamma**2 * r_diss**2))
    E_gamma_iso = corr * E_gamma_iso_star

    E_nu_pi = np.sqrt((3 * np.pi * (m_pion**5) * c**5 * Gamma**2) /
                      (8 * (m_e**2) * sigma_T * (B**2) * tau_pion)) / eV_to_erg * 1e-18
    E_nu_mu = np.sqrt((3 * np.pi * (m_muon**5) * c**5 * Gamma**2) /
                      (8 * (m_e**2) * sigma_T * (B**2) * tau_muon)) / eV_to_erg * 1e-18

    p_results = maximum_proton_energy(Gamma, r_diss, B, eps_pk_eV, alpha, beta, A, eta_acc)
    nu_results = neutrino_spectra(p_results, E_gamma_iso, xi_p, d_L_Mpc)

    print(f"\n{'='*90}\n{name} Results\n{'='*90}")
    print(f"L_iso: {L_iso:.3e},  B: {B:.3e},  Eνπ: {E_nu_pi:.5f} EeV,  Eνμ: {E_nu_mu:.5f} EeV\n")

    return {"E_nu_obs": nu_results["E_nu_obs"],
            "E2_phi_numu_earth": nu_results["E2_phi_numu_earth"],
            "B": B, "E_nu_pi": E_nu_pi, "E_nu_mu": E_nu_mu}



results = {}

for name, params in models.items():
    results[name] = run_model(name, params)



plt.figure(figsize=(8, 6))
for name, r in results.items():
    plt.loglog(r["E_nu_obs"] / 1.602e-3, r["E2_phi_numu_earth"], label=name)

plt.xlabel("Neutrino Energy $E_\\nu$ [GeV]", fontsize=13)
plt.ylabel(r"$E_\nu^2 \, \phi_{\nu_\mu}$ [erg cm$^{-2}$]", fontsize=13)
plt.xlim(1e3, 1e9)
plt.ylim(1e-8, 1e-3)
plt.legend()
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.title("Neutrino fluences for SGRB models (Kimura et al. 2017)")
plt.tight_layout()
plt.savefig("Figure1.pdf")
plt.show()
