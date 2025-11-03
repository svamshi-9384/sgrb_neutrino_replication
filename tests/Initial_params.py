import numpy as np
from scipy import integrate

def dn_depsilon(eps, A, eps_pk, alpha, beta):
    return A * np.where(eps < eps_pk,
                        (eps/eps_pk)**(-alpha),
                        (eps/eps_pk)**(-beta))

def compute_U(A, eps_pk, alpha, beta, eps_min, eps_max):
    integrand = lambda e: e * dn_depsilon(e, A, eps_pk, alpha, beta)
    result, _ = integrate.quad(integrand, eps_min, eps_max, limit=200)
    return result

def compute_L_iso_and_B():
    # Given parameters for EE-mod
    Epk_obs_keV = 1.0
    Gamma = 30.0
    alpha = 0.5
    beta  = 2.0
    L_obs = 3e48
    r_diss = 1e14
    c = 3e10
    eV_to_erg = 1.6e-12
    xi_B=0.1

    # Photon energies (units of eV)
    eps_pk       = (Epk_obs_keV / Gamma) * 1000.0
    eps_min_obs  = (0.3  / Gamma) * 1000.0
    eps_max_obs  = (10.0 / Gamma) * 1000.0
    eps_min_full = 0.1
    eps_max_full = 1e6

    # Step 1: U* in observed band
    U_obs_norm = compute_U(1.0, eps_pk, alpha, beta,
                           eps_min_obs, eps_max_obs)

    # Step 2: U over full range
    U_full_norm = compute_U(1.0, eps_pk, alpha, beta,
                            eps_min_full, eps_max_full)

    # Bolometric correction
    corr = U_full_norm / U_obs_norm

    # Normalization A
    coeff = 4 * np.pi * c * Gamma**2 * r_diss**2 * \
            U_obs_norm * eV_to_erg
    A = L_obs / coeff
    A=A/eV_to_erg

    # Compute L_iso over full range
    U_full_actual = A * U_full_norm * eV_to_erg
    L_iso = 4 * np.pi * c * Gamma**2 * r_diss**2 * U_full_actual * eV_to_erg
    B=np.sqrt(2*L_iso*xi_B/(c*((Gamma)**2)*((r_diss)**2)))
    return L_obs,corr,A,L_iso,B

def calculate_E_gamma_iso(L_obs):
    E_star_gamma_iso=1e51
    time_var=E_star_gamma_iso/L_obs  
    E_gamma_iso=L_iso*time_var   # Both these equations sum up to E_gamma_iso = Bolometriic correction * E_star_gamma_iso
    return E_gamma_iso

def calculate_E_p_iso(clf,E_gamma_iso):
    E_p_iso=clf*E_gamma_iso
    return E_p_iso

if __name__ == "__main__":
    L_obs,corr,A,L_iso,B = compute_L_iso_and_B()
    E_gamma_iso = calculate_E_gamma_iso(L_obs)
    clf=float(input("Input the value for Cosmic loading factor: "))
    E_p_iso=calculate_E_p_iso(clf,E_gamma_iso)
    
    print(f"Bolometric correction = {corr:.2f}")
    print(f"Normalization A = {A:.5e} / erg cm^-3")
    print(f"L_iso (calculated) = {L_iso:.5e} erg/s")
    print(f"B (calculated) = {B:.5e} Gauss")
    print(f"E_gamma_iso (calculated) = {E_gamma_iso:.5e} erg")
    print(f"E_p_iso (calculated) = {E_p_iso:.5e} erg")


   
