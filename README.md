# Spatially varying distributed lag non-linear models using Laplacian P-splines
Sara Rutten, Thomas Neyens, Elisa Duarte, Antonio Gasparrini and Christel Faes

## About this repository
This repository contains the Rcodes used to generate the results from the paper "Spatially varying distributed lag non-linear models using Laplacian P-splines".

## Data
Following datasets are necessary to run the simulation code and are included in the data folder:
| Dataset | Description | Downloaded from |
| --- | --- | --- |
| daily_data.RData | Daily summer temperature in Barcelona |  https://github.com/marcosqz/sbdlnm_smallarea_casestudy/tree/main/input |
| shapefile_bcn.shp | Shapefile Barcelona | https://github.com/marcosqz/sbdlnm_smallarea_casestudy/tree/main/input |

The datasets necessary to run the Sicily data application can be downloaded from https://drive.google.com/file/d/1quHnVTTpSOsYC8RnfdlLjRkSCstU-5n1/view?usp=sharing.

## Rcode
The main code data_application_Sicilia.R contains the Rcode of the data application.   
Moreover, the folder **Simulations** includes the Rcode that can be used to recover the simulation results. This folder contains several files:
| File | Description |
| --- | --- | 
| spatial_vaying_DLNM_simulation.R | Simulation code for the main simulation study |
| spatial_varying_DLNM_simulation_NB.R | Simulation code for the negative binomial simulation study |
| spatial_varying_DLNM_comparison.R | Simulation code for comparison with BAM and INLA |

The folder **functions** includes the Rcode that is needed to fit the spatially varying DLNM model. The folder contains the following files:
| File | Description |
| --- | --- | 
| DLNM_Laplace_spatially_structured.R | Function to fit the spatially varying DLNM model |
| DLNM_Laplace_spatially_structured_poisson.R | Function called by DLNM_Laplace_spatially_structured.R (fit the DLNM model with poisson distribution) |
| DLNM_Laplace_spatially_structured_NB.R | Function called by DLNM_Laplace_spatially_structured.R (fit the DLNM model with negative binomial distribution) |
| predRR_spat.R | Function to calculate the estimated RR |
| help_functions.R | Functions called by DLNM_Laplace_spatially_structured.R and predRR_spat.R |
| af_Laplace.R | Function to calculate attributable fraction |

The folder **functions Simulation** contains additional functions that are used in the simulation study. The folder contains the following files:
| File | Description |
| --- | --- | 
| crossbasis_INLA.R | Function to construct crossbasis for prediction with INLA model |
| sim_coefficients.R | Function to simulate the region-specific exposure-lag-response surface (large heterogeneity) |
| sim_coefficients_small.R | Function to simulate the region-specific exposure-lag-response surface (small heterogeneity) |
