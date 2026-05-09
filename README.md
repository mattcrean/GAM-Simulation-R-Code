# Evaluating Significance Tests in Generalised Additive Models

This repository contains the R code used for the simulation studies in my undergraduate project on the reliability of approximate significance tests for smooth terms in generalised additive models, with particular focus on the implementation in the `mgcv` package (version used at time of study was v1.8-42). 
The code investigates type I error calibration and power across several controlled scenarios, including different response distributions, sample size, smoothing parameter estimation methods, correlated covariates, and a residual permutation comparison.

All simulations were run using R v4.3.1. 
Random seeds are set before each simulation to allow for reproducibility.

## Repository structure

| File | Description |
|---|---|
| `null_sim1.R` | Type I error Simulation Scenario 1: single null smooth. Includes variation in sample size, smoothing parameter estimation method, basis dimension, and a GLM comparison. |
| `null_sim2.R` | Type I error Simulation Scenario 2: two independent null smooths. Includes sample-size comparisons, one vs two smooth comparisons, and response parameter variation. |
| `null_sim3.R` | Type I error Simulation Scenario 3: null smooth in the presence of a correlated linear covariate. |
| `null_sim4.R` | Type I error Simulation Scenario 4: two null smooths in the presence of a linear effect, with a smooth of the same covariate with linear component of spline basis removed. |
| `null_sim4_residual_permutation.R` | Residual permutation comparison for the Gaussian version of Simulation Scenario 4. |
| `power_sim1.R` | Power Simulation Scenario 1: baseline single smooth under a false null. |
| `power_sim2.R` | Power Simulation Scenario 2: nonlinear smooth effect in the presence of a linear effect of the same covariate, reflecting type I error Simulation Scenario 4. |

## Notes

The scripts are intended to reproduce the simulation results and figures discussed in the project. 
Some simulations, especially the residual permutation comparison, may take a long time to run because the GAM is refitted many times within each simulation replication.
