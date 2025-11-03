Absolutely 👍 — here’s the complete `README.md` file (ready to copy-paste directly into your repository root).

Save it as:

```bash
nano README.md
```

Paste everything below, then press `Ctrl+O`, `Enter`, `Ctrl+X` to save and exit.

---

```markdown
# SGRB Neutrino Replication (Kimura et al. 2017)

![GitHub last commit](https://img.shields.io/github/last-commit/svamshi-9384/sgrb_neutrino_replication)
![GitHub repo size](https://img.shields.io/github/repo-size/svamshi-9384/sgrb_neutrino_replication)
![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)

---

## 🧭 Overview
This repository contains an **independent computational replication** of  
**Kimura et al. (2017), _ApJ 851 – “High-Energy Neutrino Emission from Short Gamma-Ray Bursts”_**,  
focusing on the neutrino production mechanisms and spectra for five short-GRB models:

- **EE-mod** – Extended Emission (moderate)  
- **EE-opt** – Extended Emission (optimistic)  
- **Prompt**  
- **Flare**  
- **Plateau**

The calculations reproduce the spectral hierarchy and cooling-break structure discussed in the paper using a Δ-resonance approximation for photomeson production.

---

## 📂 Repository Structure
```

SGRB_neutrino_replication/
├── data/                     # Input datasets (e.g., IceCube effective area)
├── docs/                     # Section-2 summary PDF, plots
├── results/                  # Output spectra and comparison tables
├── tests/                    # Main Python scripts for each model
│   ├── sec2_neutrino_spectrum.py
│   ├── A_eff_classify.py
│   └── constants.py
├── environment.yml           # Optional conda environment
├── requirements.txt          # pip dependencies
├── LICENSE
└── README.md

````

---

## ⚙️ Installation
Clone the repository and install the required Python packages.

```bash
git clone git@github.com:svamshi-9384/sgrb_neutrino_replication.git
cd sgrb_neutrino_replication
pip install -r requirements.txt
````

Or using conda:

```bash
conda env create -f environment.yml
conda activate sgrb_neutrino
```

---

## ▶️ Usage

Run Section-2 neutrino-spectrum calculations:

```bash
python3 tests/sec2_neutrino_spectrum.py
```

Outputs:

* Fluence spectra for all five models (in `results/fluence_plots/`)
* Comparison table (`results/comparison_table.csv`)

Each run reproduces the expected hierarchy
`EE-mod > Flare > Prompt > Plateau` and
shifts (E_{\nu,\text{peak}}\propto\Gamma^2).

---

## 🧩 Next Steps

Planned extensions (2025 – 2026):

1. **Section 3:** IceCube / IceCube-Gen2 detection probabilities.
2. **GEANT4 module:** full photomeson cross-section implementation.
3. **Application to GW170817 / GRB 170817A:** off-axis vs on-axis neutrino detectability.

---

## 🧠 Physics Summary

* **Approximation:** Δ-resonance for pγ interactions
  (σ_{pγ}=5×10^{-28},\text{cm}^2,;κ_{pγ}=0.2)
* **Cooling Timescales:** (t_\text{acc}, t_\text{syn}, t_{pγ}, t_\text{dyn})
* **Main trend:** photomeson cooling dominant for EE models; adiabatic loss dominant at high Γ.
* **Peak energies:** (E_{\nu,\text{peak}}\sim10^{6–7}) GeV; fluence ≈ 10⁻⁵–10⁻⁴ erg cm⁻² at 300 Mpc.

---

## 📜 Reference

Kimura, S. S., Murase, K., Mészáros, P. (2017).
*High-Energy Neutrino Emission from Short Gamma-Ray Bursts: Prospects for Coincident Detection with Gravitational Waves.*
**ApJ 851, L55.** [https://arxiv.org/abs/1708.07075](https://arxiv.org/abs/1708.07075)

---

## 👤 Author

**Surya Vamshi Allada**
Independent Researcher, IISER Thiruvananthapuram
Email: [svamshi9384@gmail.com](mailto:svamshi9384@gmail.com)
GitHub: [@svamshi-9384](https://github.com/svamshi-9384)

---

## 🪪 License

This work is released under the **MIT License**.
Feel free to use, modify, and cite with appropriate attribution.

````

---

After saving, run:
```bash
git add README.md
git commit -m "Added professional README file"
git push
````

✅ You’ll then see the formatted README appear beautifully on your GitHub repo front page.

Would you like me to make a matching `.gitignore` next (so large data files or cache folders aren’t pushed)?
