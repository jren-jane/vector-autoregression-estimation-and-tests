# vector-autoregression
This repository carries out Matlab-coded tests for Vector Autoregression models.

## File introduction
- chowtest_chisqure_bootstrap.m: conducts Matlab-coded Chow breakpoint test and Chow sample slit test on a K-dimensional VAR(p). Those tests are used for testing structural breaks. My function returns p-values based on the asymptotic chi-square-distribution as well as bootstrap p-values based on bootstrap replications. 

- find_lag_AIC.m: I select the optimal number of lags based on AIC criterion. 

- MonteCarlo.m: I simulate 500 sets of time series data, estimate a VAR with intercept, and compute relative rejection frequencies based on the p-values obtained earlier. 

- relative_rejection_frequencies: This is the main file which returns the relative rejection frequencies.
