"""Physical constants in CGS units"""

class PhysicalConstants:
    # Masses [MeV/c^2]
    m_p_MeV = 938.272  # proton
    m_e_MeV = 0.5109989  # electron
    m_pi_MeV = 139.57  # charged pion
    m_mu_MeV = 105.658  # muon
    
    # Lifetimes [s]
    tau_pi = 2.6e-8
    tau_mu = 2.197e-6
    
    # Other constants
    c = 2.99792458e10  # cm/s
    sigma_T = 6.6524e-25  # cm^2
    e = 4.80320425e-10  # esu
    
    # Conversion factors
    MeV_to_erg = 1.60218e-6
    eV_to_erg = 1.60218e-12
    GeV_to_erg = 1.60218e-3
    Mpc_to_cm = 3.08567758e24
    keV_to_eV = 1e3
