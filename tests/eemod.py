import numpy as np
import matplotlib.pyplot as plt

c = 2.99792458e10          
m_p = 1.67262192e-24      
m_p_GeV = 0.93827208       
m_e = 9.1093897e-28        
m_pion = 2.48832e-25      
m_muon = 1.88353e-25       
e_charge = 4.80320425e-10  
sigma_T = 6.6524587e-25    
tau_pion = 2.6e-8          
tau_muon  = 2.197e-6      
eV_to_erg = 1.60217663e-12 
GeV_to_erg = 1.60217663e-3 
sigma_Delta  = 5e-28       
kappa_Delta  = 0.2        
eps_th_GeV   = 0.145       
eps_th_erg   = eps_th_GeV * GeV_to_erg

def photon_number_density(eps_prime_erg, A, eps_pk_erg, alpha, beta):

    ratio = eps_prime_erg / eps_pk_erg
    if np.isscalar(eps_prime_erg):
        if eps_prime_erg < eps_pk_erg:
            return A * ratio**(-alpha)
        else:
            return A * ratio**(-beta)

def t_acc(eps_p_erg, B_gauss, eta=1.0):
    t_acc=eta * eps_p_erg / (e_charge * c * B_gauss)
    return t_acc

def t_syn(m,eps_p_erg, B_gauss):
    t_syn = (6.0 * np.pi * (m**4) * (c**3))/((m_e**2) * sigma_T * (B_gauss**2) * eps_p_erg)
    return t_syn

def t_pgamma(eps_p_erg, eps_pk_eV, alpha, beta, A):

    eps_p_arr = np.atleast_1d(eps_p_erg).astype(float)
    gamma_p_arr = eps_p_arr / (m_p * c**2)

    eps_pk_erg = eps_pk_eV * eV_to_erg  
    epsbar_min = eps_th_erg               
    epsbar_max = 10.0 * GeV_to_erg        
    epsbar_grid = np.logspace(
        np.log10(epsbar_min),
        np.log10(epsbar_max),
        200
    )

    def photon_tail_integral_epsmin(eps0):

        eps0 = np.atleast_1d(eps0)
        I_out = np.zeros_like(eps0)

        # case1: eps0 < eps_pk
        m1 = eps0 < eps_pk_erg
        if np.any(m1):
            e0 = eps0[m1]

            I1 = (A * eps_pk_erg**(alpha) *
                  ( e0**(-(1.0+alpha)) - eps_pk_erg**(-(1.0+alpha)) ) /
                  (1.0 + alpha))
            I2 = (A * eps_pk_erg**(beta) *
                  ( eps_pk_erg**(-(1.0+beta)) ) /
                  (1.0 + beta))

            I_out[m1] = I1 + I2


        m2 = ~m1
        if np.any(m2):
            e0 = eps0[m2]
            I_high = (A * eps_pk_erg**(beta) *
                      ( e0**(-(1.0+beta)) ) /
                      (1.0 + beta))
            I_out[m2] = I_high

        return I_out

    t_pgamma_list = []

    for gamma_p in gamma_p_arr:

        eps0_grid = epsbar_grid / (2.0 * gamma_p)
        I_vals = photon_tail_integral_epsmin(eps0_grid)
        integrand = sigma_Delta * kappa_Delta * epsbar_grid * I_vals
        outer_int = np.trapz(integrand, epsbar_grid)
        t_inv = (c / (2.0 * gamma_p**2)) * outer_int
        if t_inv > 0:
            t_pgamma_list.append(1.0 / t_inv)
        else:
            t_pgamma_list.append(1e99)

    t_pgamma_arr = np.array(t_pgamma_list)
    if np.isscalar(eps_p_erg):
        return float(t_pgamma_arr[0])
    return t_pgamma_arr


def f_sup_pion(eps_p_erg, B_gauss, Gamma, r_diss):
    eps_pi_erg = 0.2 * eps_p_erg
    gamma_pi = eps_pi_erg / (m_pion * c**2)

    t_pi_dec = gamma_pi * tau_pion  

    t_pi_syn = t_syn(m_pion,eps_pi_erg, B_gauss)
    t_dyn =  r_diss / (c * Gamma)

    t_pi_cool_inv = (1.0 / t_pi_syn) + (1.0 / t_dyn)
    t_pi_cool = 1.0 / t_pi_cool_inv

    f_sup_pi = 1.0 - np.exp(-t_pi_cool / t_pi_dec)
    return f_sup_pi, t_pi_cool, t_pi_dec

def f_sup_muon(eps_p_erg, B_gauss, Gamma, r_diss):
    eps_mu_erg = 0.1 * eps_p_erg
    gamma_mu = eps_mu_erg / (m_muon * c**2)

    t_mu_dec = gamma_mu * tau_muon

    t_mu_syn = t_syn(m_muon,eps_mu_erg, B_gauss)
    t_dyn =  r_diss / (c * Gamma)

    t_mu_cool_inv = (1.0 / t_mu_syn) + (1.0 / t_dyn)
    t_mu_cool = 1.0 / t_mu_cool_inv

    f_sup_mu = 1.0 - np.exp(-t_mu_cool / t_mu_dec)
    return f_sup_mu, t_mu_cool, t_mu_dec

def f_pgamma_array(eps_p_array, t_pgamma_array, t_dyn):
    ratio = t_dyn / t_pgamma_array
    return np.minimum(1.0, ratio)


def proton_spectrum(eps_p_array, E_p_iso, eps_p_min, eps_p_max):
    norm = E_p_iso / np.log(eps_p_max / eps_p_min)
    dN_dE = norm * (eps_p_array**(-2.0))
    return dN_dE 

def neutrino_from_pion_decay(eps_p_array,
                             dNpdE,
                             f_pgamma_array,
                             f_sup_pi_array,
                             Gamma,
                             d_L_cm):

    Ep2_dN = (eps_p_array**2) * dNpdE
    E2_dNnu_src = 0.125 * f_pgamma_array * f_sup_pi_array * Ep2_dN  

    eps_nu_prime = 0.05 * eps_p_array         
    E_nu_obs = Gamma * eps_nu_prime            
    dN_dEobs = E2_dNnu_src / (E_nu_obs**2)     
    phi_nu = dN_dEobs / (4.0 * np.pi * d_L_cm**2) 

    E2_phi = (E_nu_obs**2) * phi_nu           

    return E_nu_obs, phi_nu, E2_phi  

def neutrino_from_muon_decay(eps_p_array,
                             dNpdE,
                             f_pgamma_array,
                             f_sup_pi_array,
                             f_sup_mu_array,
                             Gamma,
                             d_L_cm):

    Ep2_dN = (eps_p_array**2) * dNpdE
    E2_dNnu_src = 0.125 * f_pgamma_array * f_sup_pi_array * f_sup_mu_array * Ep2_dN

    eps_nu_prime = 0.05 * eps_p_array
    E_nu_obs = Gamma * eps_nu_prime

    dN_dEobs = E2_dNnu_src / (E_nu_obs**2)
    phi_nu = dN_dEobs / (4.0 * np.pi * d_L_cm**2)

    E2_phi = (E_nu_obs**2) * phi_nu

    return E_nu_obs, phi_nu, E2_phi  

def apply_oscillations(phi_nue_src, phi_numu_src):

    phi_e0 = phi_nue_src
    phi_mu0 = phi_numu_src
    phi_tau0 = 0.0

    phi_e_earth  = (10.0/18.0)*phi_e0 + (4.0/18.0)*(phi_mu0+phi_tau0)
    phi_mu_earth = (4.0/18.0)*phi_e0  + (7.0/18.0)*(phi_mu0+phi_tau0)
    phi_tau_earth= (4.0/18.0)*phi_e0  + (7.0/18.0)*(phi_mu0+phi_tau0)

    return phi_e_earth, phi_mu_earth, phi_tau_earth


def maximum_proton_energy(Gamma, r_diss, B, eps_pk_eV, alpha, beta, A,
                               eta_acc=1.0, n_points=800):


    eps_p_min = 10.0 * m_p * c**2
    eps_p_max_scan = 1e15 * m_p * c**2
    eps_p_array = np.logspace(np.log10(eps_p_min),
                              np.log10(eps_p_max_scan),
                              n_points)

    t_dyn =  r_diss / (c * Gamma)
    t_dyn_inv = 1.0 / t_dyn
    t_acc_array    = t_acc(eps_p_array, B, eta=eta_acc)
    t_syn_array    = t_syn(m_p,eps_p_array, B)
    t_pgamma_array = t_pgamma(eps_p_array, eps_pk_eV, alpha, beta, A)

    t_acc_inv_array    = 1.0 / t_acc_array
    t_syn_inv_array    = 1.0 / t_syn_array
    t_pgamma_inv_array = 1.0 / t_pgamma_array

    t_cool_inv_array = t_dyn_inv + t_syn_inv_array + t_pgamma_inv_array
    t_cool_array = 1.0 / t_cool_inv_array
    f_p_gamma_array = f_pgamma_array(eps_p_array, t_pgamma_array, t_dyn)
    rate_diff = t_acc_inv_array - t_cool_inv_array
    crossing_indices = np.where(np.diff(np.sign(rate_diff)))[0]

    i = crossing_indices[0]
    x1, x2 = eps_p_array[i], eps_p_array[i+1]
    y1, y2 = rate_diff[i], rate_diff[i+1]
    eps_p_max = x1 + (x2-x1)*(-y1)/(y2-y1)

    E_p_max_eV  = eps_p_max * Gamma / eV_to_erg
    E_p_max_EeV = E_p_max_eV / 1e18
    t_acc_max    = t_acc(eps_p_max, B, eta=eta_acc)
    t_syn_max    = t_syn(m_p,eps_p_max, B)
    t_pgamma_max = t_pgamma(eps_p_max, eps_pk_eV, alpha, beta, A)
    t_cool_inv_max = t_dyn_inv + (1.0/t_syn_max) + (1.0/t_pgamma_max)
    t_cool_max = 1.0 / t_cool_inv_max

    print("\n" + "="*80)
    print("RESULTS: MAXIMUM PROTON ENERGY")
    print("="*80)
    print(f"\nε_p,max' (comoving) = {eps_p_max/eV_to_erg:.3e} eV "
          f"= {eps_p_max/GeV_to_erg:.3e} GeV")
    print(f"E_p,max (observer)   = {E_p_max_EeV:.2f} EeV")


    cooling_rates = {
        'Adiabatic': t_dyn_inv,
        'Synchrotron': 1.0/t_syn_max,
        'Photomeson': 1.0/t_pgamma_max
    }
    dominant = max(cooling_rates, key=cooling_rates.get)
    print(f"\nDominant cooling mechanism at Ep,max: {dominant}")
    for name, rate in cooling_rates.items():
        pct = 100.0 * rate / t_cool_inv_max
        print(f"  {name:12s}: {rate:.3e} s^-1 ({pct:5.1f}%)")

    return {
        'Gamma': Gamma,
        'r_diss': r_diss,
        'B': B,
        'alpha': alpha,
        'beta': beta,
        'A': A,
        'eps_pk_eV': eps_pk_eV,
        'eps_p_array': eps_p_array,
        'eps_p_min': eps_p_min,
        'eps_p_max': eps_p_max,
        't_dyn': t_dyn,
        't_dyn_inv': t_dyn_inv,
        't_acc_array': t_acc_array,
        't_syn_array': t_syn_array,
        't_pgamma_array': t_pgamma_array,
        't_acc_inv_array': t_acc_inv_array,
        't_syn_inv_array': t_syn_inv_array,
        't_pgamma_inv_array': t_pgamma_inv_array,
        't_cool_inv_array': t_cool_inv_array,
        't_cool_array': t_cool_array,
        'f_p_gamma_array': f_p_gamma_array,
        'E_p_max_EeV': E_p_max_EeV
    }

def neutrino_spectra(proton_results,
                               E_gamma_iso,
                               xi_p,
                               d_L_Mpc):



    eps_p_array    = proton_results['eps_p_array']
    eps_p_min      = proton_results['eps_p_min']
    eps_p_max      = proton_results['eps_p_max']
    f_p_gamma_array= proton_results['f_p_gamma_array']
    Gamma          = proton_results['Gamma']
    r_diss         = proton_results['r_diss']
    B              = proton_results['B']
    t_dyn      = proton_results['t_dyn']

    E_p_iso = xi_p * E_gamma_iso
    d_L_cm = d_L_Mpc * 3.086e24  # Mpc -> cm

    dN_p_dE_p = proton_spectrum(eps_p_array, E_p_iso, eps_p_min, eps_p_max)

    f_sup_pi_array = np.zeros_like(eps_p_array)
    f_sup_mu_array = np.zeros_like(eps_p_array)

    for i, Ep_prime in enumerate(eps_p_array):
        f_sup_pi, _, _ = f_sup_pion(Ep_prime, B, Gamma, r_diss)
        f_sup_mu, _, _ = f_sup_muon(Ep_prime, B, Gamma, r_diss)
        f_sup_pi_array[i] = f_sup_pi
        f_sup_mu_array[i] = f_sup_mu
    

    E_nu_obs_pi, phi_numu_src_pi, E2_phi_numu_src_pi = neutrino_from_pion_decay(
        eps_p_array,
        dN_p_dE_p,
        f_p_gamma_array,
        f_sup_pi_array,
        Gamma,
        d_L_cm
    )


    E_nu_obs_mu, phi_nu_from_mu_src, E2_phi_nu_from_mu_src = neutrino_from_muon_decay(
        eps_p_array,
        dN_p_dE_p,
        f_p_gamma_array,
        f_sup_pi_array,
        f_sup_mu_array,
        Gamma,
        d_L_cm
    )

    phi_nue_src  = phi_nu_from_mu_src
    phi_numu_src = phi_numu_src_pi + phi_nu_from_mu_src


    phi_nue_earth, phi_numu_earth, phi_nutau_earth = apply_oscillations(
        phi_nue_src,
        phi_numu_src
    )
    

    E2_phi_nue_earth   = (E_nu_obs_pi**2) * phi_nue_earth
    E2_phi_numu_earth  = (E_nu_obs_pi**2) * phi_numu_earth
    E2_phi_nutau_earth = (E_nu_obs_pi**2) * phi_nutau_earth

    idx_peak_mu = np.argmax(E2_phi_numu_earth)
    idx_peak_e  = np.argmax(E2_phi_nue_earth)
    Ep_to_Enu_factor = 0.05 * Gamma / eV_to_erg / 1e18  
    print("\n" + "="*80)
    print("NEUTRINO RESULTS (after oscillations)")
    print("="*80)

    print(f"\nPeak ν_mu (Earth) E^2 φ at:")
    print(f"  E_ν,peak ~ {E_nu_obs_pi[idx_peak_mu]/eV_to_erg/1e9:.3e} GeV")
    print(f"  E^2 φ_νμ ~ {E2_phi_numu_earth[idx_peak_mu]:.3e} erg/cm^2")

    print(f"\nPeak ν_e (Earth) E^2 φ at:")
    print(f"  E_ν,peak ~ {E_nu_obs_pi[idx_peak_e]/eV_to_erg/1e9:.3e} GeV")
    print(f"  E^2 φ_νe ~ {E2_phi_nue_earth[idx_peak_e]:.3e} erg/cm^2")

    total_phi = (phi_nue_earth[idx_peak_mu] +
                 phi_numu_earth[idx_peak_mu] +
                 phi_nutau_earth[idx_peak_mu])
    fe = phi_nue_earth[idx_peak_mu]   / total_phi if total_phi>0 else 0
    fmu= phi_numu_earth[idx_peak_mu]  / total_phi if total_phi>0 else 0
    ftau=phi_nutau_earth[idx_peak_mu] / total_phi if total_phi>0 else 0



    return {
        'E_nu_obs': E_nu_obs_pi,  # erg
        'phi_nue_earth': phi_nue_earth,
        'phi_numu_earth': phi_numu_earth,
        'phi_nutau_earth': phi_nutau_earth,
        'E2_phi_nue_earth': E2_phi_nue_earth,
        'E2_phi_numu_earth': E2_phi_numu_earth,
        'E2_phi_nutau_earth': E2_phi_nutau_earth,
        'f_sup_pi_array': f_sup_pi_array,
        'f_sup_mu_array': f_sup_mu_array,
        'f_p_gamma_array': f_p_gamma_array,
        'dN_p_dE_p': dN_p_dE_p,
        'eps_p_array': eps_p_array,
        'Gamma': Gamma
    }
    



def plot_timescales_and_fp(results):

    eps_p_array        = results['eps_p_array']
    t_acc_inv_array    = results['t_acc_inv_array']
    t_syn_inv_array    = results['t_syn_inv_array']
    t_pgamma_inv_array = results['t_pgamma_inv_array']
    t_cool_inv_array   = results['t_cool_inv_array']
    t_dyn_inv          = results['t_dyn_inv']
    f_p_gamma_array    = results['f_p_gamma_array']
    Gamma              = results['Gamma']

    E_p_obs_eV  = eps_p_array * Gamma / eV_to_erg
    E_p_obs_EeV = E_p_obs_eV / 1e18

    fig, axes = plt.subplots(2, 1, figsize=(12, 10))


    ax = axes[0]
    ax.loglog(E_p_obs_EeV, t_acc_inv_array, linewidth=2.0, label=r'$t_{\mathrm{acc}}^{-1}$')
    ax.loglog(E_p_obs_EeV, t_cool_inv_array, linewidth=2.0, label=r'$t_{\mathrm{cool}}^{-1}$')
    ax.loglog(E_p_obs_EeV, t_syn_inv_array, linestyle='--', linewidth=1.5, label=r'$t_{\mathrm{syn}}^{-1}$')
    ax.loglog(E_p_obs_EeV, t_pgamma_inv_array, linestyle='--', linewidth=1.5, label=r'$t_{p\gamma}^{-1}$')
    ax.hlines(t_dyn_inv, E_p_obs_EeV.min(), E_p_obs_EeV.max(),
              linestyles=':', linewidth=1.5, label=r'$t_{\mathrm{dyn}}^{-1}$')
    ax.set_xlabel(r'$E_p \, (\mathrm{EeV})$')
    ax.set_ylabel(r'Rate [s$^{-1}$]')
    ax.set_title("Acceleration vs Cooling Rates")
    ax.grid(True, which='both', alpha=0.3)
    ax.legend(loc='best')

    ax2 = axes[1]
    ax2.semilogx(E_p_obs_EeV, f_p_gamma_array, linewidth=2.0, label=r'$f_{p\gamma}$')
    ax2.set_xlabel(r'$E_p \, (\mathrm{EeV})$')
    ax2.set_ylabel(r'$f_{p\gamma}$')
    ax2.set_ylim(0, 1.05)
    ax2.set_title("Meson Production Efficiency")
    ax2.grid(True, which='both', alpha=0.3)
    ax2.legend(loc='best')

    plt.tight_layout()
    plt.show()

def plot_neutrino_spectra(neu_results):

    E_nu_obs = neu_results['E_nu_obs']              # erg
    E2_phi_nue   = neu_results['E2_phi_nue_earth']   # erg/cm^2
    E2_phi_numu  = neu_results['E2_phi_numu_earth']  # erg/cm^2
    E2_phi_nutau = neu_results['E2_phi_nutau_earth'] # erg/cm^2

    E_nu_GeV = E_nu_obs / eV_to_erg / 1e9

    plt.figure(figsize=(8,6))
    plt.loglog(E_nu_GeV, E2_phi_nue,   label=r'$\nu_e$ (Earth)')
    plt.loglog(E_nu_GeV, E2_phi_numu,  label=r'$\nu_\mu$ (Earth)')
    plt.loglog(E_nu_GeV, E2_phi_nutau, label=r'$\nu_\tau$ (Earth)')

    plt.xlabel(r'$E_\nu \, [\mathrm{GeV}]$')
    plt.ylabel(r'$E_\nu^2 \Phi_\nu \, [\mathrm{erg} \, \mathrm{cm}^{-2}]$')
    plt.title("Predicted Neutrino Fluence at Earth")
    plt.grid(True, which='both', alpha=0.3)
    plt.legend(loc='best')
    plt.tight_layout()
    plt.show()



if __name__ == "__main__":

    Gamma        = 30.0
    r_diss       = 1.0e14          
    B            = 2.920420e3      
    alpha        = 0.5
    beta         = 2.0


    eps_pk_comoving_erg = 5.341e-11  
    eps_pk_eV = eps_pk_comoving_erg / eV_to_erg  
    A = 1.73920e+14/eV_to_erg  # [photons / (erg cm^3)]
    # Burst energetics etc.
    L_gamma_iso = 1.200e49    
    E_gamma_iso = 3.800e51    
    xi_p        = 10.0        
    d_L_Mpc     = 300.0        


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

    neu_results = neutrino_spectra(
        proton_results,
        E_gamma_iso=E_gamma_iso,
        xi_p=xi_p,
        d_L_Mpc=d_L_Mpc
    )

    # Step 3: plots
    plot_timescales_and_fp(proton_results)
    plot_neutrino_spectra(neu_results)




