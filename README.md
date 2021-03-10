# vector-autoregression
This repository carries out Matlab-coded estimation and tests for Vector Autoregression models.

## File introduction
- chowtest_chisqure_bootstrap.m: This is a function that conducts Chow breakpoint test and Chow sample slit test on a K-dimensional VAR(p). Those tests are used for testing structural breaks. It takes as inputs a (T+p)K matrix, a break point, the number of lags, and an indicator of whether to include an intercept. It returns p-values based on the asymptotic chi-square-distribution as well as bootstrap p-values based on bootstrap replications. 
  - In the first half of the function, the model is estimated, the values of the BP and SS test statistics are obtained and stored, and their $p$-values are computed using an asymptotic chi-square-distribution. 
  - In the second half of the function, bootstrapping is implemented, values of the test statistics from bootstrapping are collected, and the original test statistics from the first half are compared with the bootstrapped distribution to show the p-values.

- find_lag_AIC.m: This is a function that takes as inputs the (T+p)K matrix, the maximum number of lags, and an indicator of whether to include an intercept, and gives as outputs the lag order chosen by the AIC information criterion.

- MonteCarlo.m: It takes as inputs the number of simulations, the number of observations to be discarded, and the length of the time series. It gives as outputs all the replicated time series in one matrix.  

- relative_rejection_frequencies: This is my main file which specifies the maximum number of lags, the number of simulations, the number of observations to be discarded, an indicator of including an intercept, and the significance level. It then specifies the different choices for the length of data and the break point. 
