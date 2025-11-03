#!/usr/bin/env python3
"""
make_Aeff_auto.py
------------------------------------------------------------
Reads IceCube IC86-I effective area table and automatically computes:
  1️⃣ Aeff_up.txt   (declination < 0°)
  2️⃣ Aeff_down.txt (declination ≥ 0°)
  3️⃣ Aeff_all.txt  (declination-averaged, full sky)

Usage:
  Just run:
      python make_Aeff_auto.py
------------------------------------------------------------
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# =============================================================================
# STEP 1: Load the input data
# =============================================================================
def load_ic86_data(path="IC86_I_effectiveArea.csv"):
    try:
        df = pd.read_csv(
            path,
            delim_whitespace=True,
            comment="#",
            names=["logEmin", "logEmax", "Decmin", "Decmax", "Aeff_cm2"],
            header=None,
            skiprows=1
        )
    except Exception as e:
        raise RuntimeError(f"❌ Could not read {path}: {e}")

    if df.empty:
        raise RuntimeError("❌ The file seems empty or malformed.")
    
    df["E_cen_GeV"] = 10 ** (0.5 * (df["logEmin"] + df["logEmax"]))
    df["Dec_cen"] = 0.5 * (df["Decmin"] + df["Decmax"])
    print(f"✅ Loaded {len(df)} rows from {path}")
    return df


# =============================================================================
# STEP 2: Hemisphere averages (up and down)
# =============================================================================
def compute_Aeff_hemisphere(df):
    print("\n🌎 Computing hemisphere-based effective areas...")
    mask_up = df["Dec_cen"] < 0
    mask_dn = df["Dec_cen"] >= 0
    E_unique = np.sort(df["E_cen_GeV"].unique())
    A_up, A_dn = [], []

    for e in E_unique:
        mE = np.isclose(df["E_cen_GeV"], e, rtol=0.05)
        up_vals = df.loc[mE & mask_up, "Aeff_cm2"].values
        dn_vals = df.loc[mE & mask_dn, "Aeff_cm2"].values
        A_up.append(np.nanmean(up_vals) if len(up_vals) else np.nan)
        A_dn.append(np.nanmean(dn_vals) if len(dn_vals) else np.nan)

    mask_valid_up = ~np.isnan(A_up)
    mask_valid_dn = ~np.isnan(A_dn)

    np.savetxt("Aeff_up.txt",
               np.c_[E_unique[mask_valid_up], np.array(A_up)[mask_valid_up]],
               fmt="%.6e", header="E_GeV  Aeff_up_cm2")
    np.savetxt("Aeff_down.txt",
               np.c_[E_unique[mask_valid_dn], np.array(A_dn)[mask_valid_dn]],
               fmt="%.6e", header="E_GeV  Aeff_down_cm2")

    print(f"✅ Wrote Aeff_up.txt   ({mask_valid_up.sum()} bins)")
    print(f"✅ Wrote Aeff_down.txt ({mask_valid_dn.sum()} bins)")
    return E_unique, np.array(A_up), np.array(A_dn), mask_valid_up, mask_valid_dn


# =============================================================================
# STEP 3: Declination-averaged (full-sky)
# =============================================================================
def compute_Aeff_allsky(df):
    print("\n🌐 Computing declination-averaged effective area...")
    mask_nonzero = df["Aeff_cm2"] > 0
    E_unique = np.sort(df.loc[mask_nonzero, "E_cen_GeV"].unique())
    Aeff_mean = []
    for e in E_unique:
        group = df[np.isclose(df["E_cen_GeV"], e, rtol=0.05) & mask_nonzero]
        Aeff_mean.append(group["Aeff_cm2"].mean())

    np.savetxt("Aeff_all.txt", np.c_[E_unique, Aeff_mean],
               fmt="%.6e", header="E_GeV  Aeff_cm2 (declination-averaged)")
    print(f"✅ Wrote Aeff_all.txt ({len(E_unique)} bins)")
    return E_unique, np.array(Aeff_mean)


# =============================================================================
# STEP 4: Plot results
# =============================================================================
def plot_results(E_up, A_up, mask_up, E_dn, A_dn, mask_dn, E_all, A_all):
    print("\n📈 Plotting results...")
    plt.figure(figsize=(7, 5))
    plt.loglog(E_up[mask_up],  A_up[mask_up]/1e4,  'b.-', label="Up-going (m²)")
    plt.loglog(E_dn[mask_dn],  A_dn[mask_dn]/1e4,  'r.-', label="Down-going (m²)")
    plt.loglog(E_all,          A_all/1e4,          'g--', label="All-sky avg (m²)")
    plt.xlabel("Neutrino energy [GeV]", fontsize=12)
    plt.ylabel("Effective area [m²]", fontsize=12)
    plt.title("IceCube IC86-I Effective Area (Hemisphere + All-sky)", fontsize=13)
    plt.grid(True, which="both", alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.show()
    print("✅ Plot complete!\n")


# =============================================================================
# MAIN
# =============================================================================
def main():
    print("🚀 Starting IceCube A_eff computation...")
    df = load_ic86_data()

    # Hemisphere averages
    E_unique, A_up, A_dn, mask_up, mask_dn = compute_Aeff_hemisphere(df)

    # All-sky average
    E_all, A_all = compute_Aeff_allsky(df)

    # Plot
    plot_results(E_unique, A_up, mask_up, E_unique, A_dn, mask_dn, E_all, A_all)

    print("🎯 All tasks done! Files ready for use in neutrino detection module (sec3.py).")


if __name__ == "__main__":
    main()
