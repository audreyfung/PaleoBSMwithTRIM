# Track Length Distribution Calculator

Minerals deep underground can retain traces of nuclear recoil damage tracks within their crystal lattice. This code calculates the **track length distribution** caused by nuclear recoils from various sources, including:

- Neutrinos  
- Neutrons from radioactive decays  
- Dark matter particles  
- Other rare sources

For a full description of the physics and methodology, please refer to:  
📄 [arXiv:2504.08885](https://arxiv.org/abs/2504.08885)

---

## Overview

This code combines results from **TRIM** (Transport of Ions in Matter) simulations — datasets available at [Zenodo](https://zenodo.org/records/15358408) — with analytical modeling of incoming particle fluxes, calibrated by observational data. The result is a prediction of the differential recoil track rate as a function of track length:

$$
\frac{dR}{dx}(x) = \sum_i \int dE_R \, P_i(x \mid E_R) \, \frac{dR_i}{dE_R}(E_R) \, \mathcal{P}_{i, \mathrm{track}}(E_R)
$$

---

## Code Structure

- `compute_trackSpectra.py`  
  Computes the track length distributions using recoil spectra $\frac{dR_i}{dE_R}(E_R)$ from various sources.

- `wimps/compute_wimps.py`  
  Computes recoil spectra from WIMP interactions.

- `light_mediators/compute_neutrino_spectra.py`  
  Computes recoil spectra from neutrino interactions mediated by light new particles.

---

## Mineral Object

Minerals are defined as Python objects of class `Minerals`.

Each `Minerals` instance includes:

- Chemical composition  
- Atomic numbers and weights  
- Atomic fractions of constituent elements  
- Path to associated SRIM simulation data

Currently, only **Olivine** is supported.

## Example Output

A sample analysis is provided in the Jupyter notebook:  
- `Example.ipynb`

This notebook illustrates example results for track length distributions from different recoil sources.

---
