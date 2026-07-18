# DFT Study of a Distorted CuO₂ Lattice

## Overview

This repository contains the input files and post-processing scripts used for Density Functional Theory (DFT) and phonon calculations on a distorted CuO₂ lattice using Quantum ESPRESSO.

The project was a part of my bachelor graduation thesis  at Cairo University.

## Objectives

- Perform self-consistent DFT calculations
- Calculate phonon properties
- Extract wavefunction (`psi`) and `dvscf` quantities
- Analyze the electronic and vibrational properties of the distorted lattice
- calculate the electron phonon coupling matrix element for distorted and undistorted lattices 
- comparing the results and see if it verifies theoretical predictions 

## Software

- Quantum ESPRESSO
- Python 3
- NumPy
- Matplotlib

## Repository Structure

```
input_files/
    scf.in
    ph.in
    q2r.in
    matdyn.in
    wfc.in
scripts/
g calc final script.py
 dvscf_txt.py
 polarization_txt.py
figures/

docs/
    workflow.md
```

# Workflow

1. Run `pw.x` using `scf.in`
2. Run `ph.x` using `ph.in`
3. Convert dynamical matrices with `q2r.x`
4. Compute phonon dispersion using `matdyn.x`
5. Run the Python scripts for post-processing

 #Future Work

- Electron–phonon coupling calculations
- Band structure analysis
- Density of states
- Improved visualization tools

# Author

**Mohamed Hossam Abdelghafar**

B.Sc. Physics

Cairo University
