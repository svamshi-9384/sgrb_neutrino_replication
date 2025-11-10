# SGRB Neutrino Replication (Kimura et al. 2017)



## 🧭 Overview
This repository contains an **independent computational replication** of  
**Kimura, Murase & Mészáros (2017), _ApJ 851 – “High-Energy Neutrino Emission from Short Gamma-Ray Bursts”_**,  
focusing on neutrino production mechanisms and fluence spectra for five short-GRB models:

- **EE-mod** – Extended Emission (moderate)  
- **EE-opt** – Extended Emission (optimistic)  
- **Prompt**  
- **Flare**  
- **Plateau**

The calculations reproduce the spectral hierarchy and trends discussed in the paper using a Δ-resonance approximation for photomeson production.

---

## 📂 Repository Structure
```

sgrb_neutrino_reproduction/
├── README.md
├── environment.yml
├── requirements.txt
│
├── data/
│   └── table1_parameters.json
│
├── docs/
│   ├── Kimura_Murase_2017.pdf
│   ├── Murase_2006.pdf
│   └── Murase_2018.pdf
│
├── notebooks/
│
├── output/
│   ├── figures/
│   ├── logs/
│   └── tables/
│
├── src/
│   ├── **init**.py
│   ├── **pycache**/
│   └── constants.py
│
└── tests/
├── Aeff.py
├── all_models.py
├── eemod.py
├── sec3.py
├── Initial_params.py
├── IC86_I_effectiveArea.csv
├── Aeff_all.txt
├── Aeff_down.txt
├── Aeff_up.txt
├── all.txt
├── extract_EE_mod_spectrum.py
└── **pycache**/

````

---

## ⚙️ Installation
Clone the repository and install dependencies.

```bash
git clone git@github.com:svamshi-9384/sgrb_neutrino_replication.git
cd sgrb_neutrino_replication
pip install -r requirements.txt
````

or with Conda:

```bash
conda env create -f environment.yml
conda activate sgrb_neutrino
```

---

## ▶️ Usage

### 1️⃣ Section 2 — Neutrino Fluence Spectrum

Run:

python3 tests/all_models.py


or, for a single model:


python3 tests/eemod.py



### 2️⃣ Section 3 — IceCube/IceCube-Gen2 Detection Probabilities

Run:


python3 tests/sec3.py


This estimates the detection probability
( P = 1 - e^{-N_\mu} ),
where
( N_\mu = \int \phi_\nu(E_\nu) A_{\text{eff}}(E_\nu) , dE_\nu )
using effective areas from `IC86_I_effectiveArea.csv`.

### Application to GRB170817A

Change parameters in Initial_params.py, then substitute the outputs to eemod.py and sec3.py.

## 🧠 Physics Summary

* **Photomeson interaction:** Δ-resonance approximation
  (σ_{pγ}=5×10^{-28},\text{cm}^2,;κ_{pγ}=0.2)
* **Cooling timescales:** (t_\text{acc}, t_\text{syn}, t_{pγ}, t_\text{dyn})
* **Spectral trends:**

  * EE models yield highest (E^2\phi_\nu)
  * (E_{\nu,\text{peak}}\propto\Gamma^2)
  * Adiabatic losses dominate at high Lorentz factors

---

## 🧩 Future Work

* Extend to **full GEANT4 photomeson cross-section**.
* Apply framework to **GW170817 / GRB 170817A** for off-axis detectability estimates.
* Publish extended results as a **full AAS paper** (2026 target).

---

## 📜 Reference

Kimura, S. S., Murase, K., & Mészáros, P. (2017).
*High-Energy Neutrino Emission from Short Gamma-Ray Bursts: Prospects for Coincident Detection with Gravitational Waves.*
**ApJ 851, L55.** [https://arxiv.org/abs/1708.07075](https://arxiv.org/abs/1708.07075)

---

## 👤 Author

**Surya Vamshi Allada**
Independent Researcher, IISER Thiruvananthapuram
Email: [svamshi9384@gmail.com](mailto:svamshi9384@gmail.com)
GitHub: [@svamshi-9384](https://github.com/svamshi-9384)

---


