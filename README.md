# vector-autoregression
This repository carries out Matlab-coded tests for Vector Autoregression models.

The chowtest_chisqure_bootstrap file conducts Matlab-coded Chow breakpoint test and Chow sample slit test on a K-dimensional VAR(p). Those breakpoint tests are used for testing structural breaks. My function returns p-values based on the asymptotic chi-square-distribution as well as bootstrap p-values based on bootstrap replications. 

In the find_lag file, I select the optimal number of lags based on AIC. 

In the MonteCarlo file, I simulate 500 sets of time series data, estimate a VAR with intercept, and compute relative rejection frequencies based on the p-values obtained earlier. 
