# Scripts accompanying Bhat and Guttal 2024

This repository contains scripts and data accompanying "Eco-evolutionary dynamics for finite populations from first principles and noise-induced reversal of selection", to appear in *The American Naturalist*. The scripts and data files in this repository can be used to reproduce the plots presented in Box 3 and Supplementary section S12 of the manuscript. A citable DOI for this repo can be found here: [![DOI](https://zenodo.org/badge/848828414.svg)](https://zenodo.org/doi/10.5281/zenodo.13384639).

## Guide to the repo

This repo contains two folders, ```box_3_example``` and ```supp_LV_example```, which contain scripts for reproducing the figures presented in Box 3 and Supplementary section S12 respectively.

### Scripts for reproducing the figure in box 3

```direct_mechanism.R``` and ```indirect_mechanism.R``` are R scripts that reproduce Box Figure IIA and Box Figure IIB in the main text respectively. The plots generated are stored as SVG files in the ```plots``` folder. Both scripts plot the speed density of the example introduced in Box 3 for various parameter values.


### Scripts for reproducing the figure in the Supplementary

The ```data``` folder contains CSV files that are output by the scripts presented in the main repository. The ```plots``` folder contains plots that make up figure S2 in the manuscript. The plots were arranged into a single figure using Inkscape. The scripts serve the following functions:

- ```noise_induced_selection_model_SSA.R``` is an R script that uses Gillespie's stochastic simulation algorithm (SSA) to simulate the exact birth-death process corresponding to the resource competition model we introduce as an example in supplementary section S12. The output of this file is a CSV file stored in the ```data``` folder.
- ```analytical_solution_in_deterministic_limit.ipynb``` is a Jupyter Notebook in which we use Python code to symbolically find the fixed points of the cubic ODE corresponding to the infinite population (deterministic) limit of our birth-death process.
- ```compare_SSA_with_numerical_integration.R``` is an R script that compares the steady state density obtained from the SSA (using ```noise_induced_selection_model_SSA.R```) with the quasi-steady state density obtained by numerically solving the Fokker-Planck equation corresponding to the Ito SDE obtained after a system-size expansion of the exact birth-death process by plotting the two density functions on the same graph. We also plot the infinite population limit prediction obtained in  ```analytical_solution_in_deterministic_limit.ipynb```. The output of this file is an SVG file stored in the ```plots``` folder. These plots are figure S2A and figure S2B in the manuscript.
