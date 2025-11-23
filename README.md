# simInter
Data and codes for "Efficient Interaction Analysis in RCTs"

## Overview

This repository contains simulation code for evaluating different randomization methods in Randomized Controlled Trials (RCTs) when analyzing treatment-covariate interactions. The project compares the performance of various randomization schemes and statistical tests for interaction effects.

## Randomization Methods

The project implements and compares the following randomization schemes:

- **SRS**: Simple Random Sampling
- **SBR**: Stratified Block Randomization
- **SBCD**: Stratified Efron's Biased Coin Design
- **PS**: Pocock and Simon's Minimization 

## Statistical Tests

The code evaluates four different test procedures:

1. **OLS Test**: Ordinary Least Squares test
2. **Usual Test**: Standard interaction test
3. **Modified Test**: Modified test accounting for randomization
4. **Efficient Test**: Efficient test using cross-fitting

## File Structure

- `main.R`: Main simulation script that orchestrates Monte Carlo simulations and generates results tables
- `randomization_methods.R`: Implementation of randomization schemes (SRS, SBR, SBCD, PS)
- `sample_model.R`: Data generation functions and model fitting procedures for different outcome models
- `test_formula.R`: Test statistic computation and output functions

## Dependencies

Required R packages:

```r
library(caratINT)      # For interaction analysis functions
library(dplyr)         # Data manipulation
library(sandwich)      # Robust standard errors
library(ranger)        # Random forest models
library(np)            # Nonparametric regression
library(doSNOW)        # Parallel computing
library(parallel)      # Parallel computing utilities
library(xtable)        # Table generation
library(bigsplines)    # Spline functions
library(caret)         # Classification and regression training
```

## Installation

1. Install required R packages:
```r
install.packages(c("dplyr", "sandwich", "ranger", "np", "doSNOW", 
                   "parallel", "xtable", "bigsplines", "caret"))
```

2. Install the `caratINT` package:
```r
if(!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
devtools::install_github('zhanglk6/caratINT')
```

## Usage

### Basic Simulation

Run a simulation for a specific model:

```r
source("main.R")

# Simulate Model 1 with N=800, delta=0, pi=0.5
table1_no <- simulate_and_compute(1, N=800, delta=0, pi=0.5, transformed=FALSE)
```

### Generate Paper Tables

Generate complete results tables for the paper:

```r
# Standard output with pi=0.5, N=800
tab_1 <- paper_output(pi=1/2, N=800)

# Output with pi=2/3
tab_2 <- paper_output(pi=2/3, N=800)

# Small sample sizes
tab_small <- paper_output_small(pi=1/2)

# Transformed covariates
tab_trans <- paper_output_trans(pi=1/2, N=800)
```

### Custom Simulations

Run custom simulations:

```r
# Simulate specific model with custom parameters
result <- simulate_model(
  model_num = 1,      # Model number (1-4)
  N = 800,            # Sample size
  times = 5000,       # Number of Monte Carlo replications
  delta = 0,          # Interaction effect size
  pi = 0.5,           # Treatment assignment probability
  seeds = 439,        # Random seed
  transformed = FALSE # Whether to transform covariates
)

# Compute rejection proportions
reject_prop <- compute_reject_prop(result, threshold = 0.05)
```

## Models

The code implements four different outcome models:

1. **Model 1**: Linear model with continuous covariates
2. **Model 2**: Exponential model
3. **Model 3**: Log-linear model
4. **Model 4**: Log-linear model with high-dimensional covariates

## Parallel Computing

The simulation uses parallel computing with 8 cores by default. Adjust the number of cores in `main.R`:

```r
numCores <- 8  # Change to desired number of cores
```

## Output

Results are returned as matrices with rejection proportions (multiplied by 100) for each combination of:
- Randomization method (SRS, SBR, SBCD, PS)
- Test procedure (OLS, Usual, Modified, Efficient)


