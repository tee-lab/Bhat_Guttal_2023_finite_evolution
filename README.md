# Scripts accompanying Bhat and Guttal 2023

This repository contains scripts and data accompanying Bhat and Guttal 2023, "Eco-evolutionary dynamics for finite populations from first principles and noise-induced reversal of selection". The scripts and data files in this repository are used to generate the plots presented in figure 2 in the manuscript.

## Guide to the repo

The ```data``` folder contains CSV files that are output by the scripts presented in the main repository. The ```plots``` folder contains plots that make up figure 2 in the manuscript. These plots were arranged into a single figure using inkscape. The scripts are to be used as follows:

- ```noise_induced_selection_model_SSA.R``` is an R script that uses Gillespie's stochastic simulation algorithm (SSA) to simulate the exact birth-death process corresponding to the example we introduce in section 4 of the manuscript. The output of this file is a CSV file stored in the ```data``` folder.
- ```analytical_solution_in_deterministic_limit.ipynb``` is a Jupyter Notebook in which we symbolically solve the cubic ODE corresponding to the infinite population (deterministic) limit of our birth-death process.
- ```compare_SSA_with_numerical_integration.R``` compares the steady state density obtained from the SSA (using ```noise_induced_selection_model_SSA.R```) with the steady state density obtained by numerically solving the Fokker-Planck equation corresponding to the Ito SDE obtained after a system-size expansion of the exact birth-death process by plotting the two density functions on the same graph. We also plot the infinite population limit prediction as obtained in  ```analytical_solution_in_deterministic_limit.ipynb```. The output of this file is an SVG file stored in the ```plots``` folder. These plots are used to make figure 2A and figure 2B in the manuscript.
- ```numerical_integration_for_param_sweep.py``` numerically solves the Fokker-Planck equation corresponding to the Ito SDE obtained after a system-size expansion of the exact birth-death process for a range of values of K and epsilon. The output of this file is a CSV file titled ```numerical_eps_K_param_sweep.csv``` stored in the ```data``` folder, as well as a heatmap of this data (as a PNG file) stored in the ```plots``` folder. This plot is figure 2C in the manuscript.