# vector-autoregression
These are Matlab-coded Chow breakpoint test and Chow sample slit test on a K-dimensional Vector Autoregression model. Those breakpoint tests are used for testing structural breaks.

My function returns p-values based on the asymptotic chi-square-distribution as well as bootstrap p-values based on bootstrap replications. For the Monte Carlo experiment, I simulate 500 sets of time series data, estimate a VAR with intercept, and compute relative rejection frequencies based on the p-values obtained earlier. 
