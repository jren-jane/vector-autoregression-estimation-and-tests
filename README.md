# vector-autoregression-estimation-and-tests
This repository carries out Matlab-coded estimation and tests for Vector Autoregression models. I explore the finite sample properties of the Chow breakpoint tests and the Chow sample split using some Monte Carlo experiments.

## File introduction
- chowtest_chisquare_bootstrap.m: This is a function that conducts Chow breakpoint test and Chow sample slit test on a K-dimensional VAR(p). Those tests are used for testing structural breaks. It takes as inputs a (T+p)K matrix, a break point, the number of lags, and an indicator of whether to include an intercept. It returns p-values based on the asymptotic chi-square-distribution as well as bootstrap p-values based on bootstrap replications. 
  - In the first half of the function, the model is estimated, the values of the BP and SS test statistics are obtained and stored, and their p-values are computed using an asymptotic chi-square-distribution. 
  - In the second half of the function, instead of assuming that the test statistic follows asymptotically a chi-square distribution, we bootstrap historical values. Values of the test statistics are collected from each bootstrap and forms our sample distribution of statistic. The original test statistics from the first half are compared with the bootstrapped distribution to show the p-values.

- find_lag_AIC.m: This is a function that takes as inputs the (T+p)K matrix, the maximum number of lags, and an indicator of whether to include an intercept, and gives as outputs the lag order chosen by the AIC information criterion.

- MonteCarlo.m: This file contains our data generating process. It takes as inputs the number of simulations, the number of observations to be discarded, and the length of the time series. It gives as outputs all the replicated time series in one matrix.  

- relative_rejection_frequencies: This is my main file which specifies the maximum number of lags, the number of simulations, the number of observations to be discarded, an indicator of including an intercept, and the significance level. It then specifies the different choices for the length of data and the break point. 

## Results (You results may vary, which depend on the actual simulation)
Case | T=80 Tb=0.5T | T=500 Tb=0.5T | T=80 Tb=0.2T | T=500 Tb=0.2T
---- | ------------ | ------------- | ------------ | -------------
BP, asymptotic chi-square | 0.15 | 0.05 | 0.24 | 0.06
SS, asymptotic chi-square | 0.14 | 0.06 | 0.16 | 0.06
BP, bootstrapped | 0.09 | 0.05 | 0.10 | 0.05
SS, bootstrapped | 0.11 | 0.06 | 0.11 | 0.04

## Implications
At a significance level of 5%, we should be expecting our relative rejection frequencies to be around 0.05. But we see many exceptions. The result shows that
- Asymptotic chi-square distribution is a poor approximation of the actual distribution.
- Shorter samples tend to give inaccurate test results. 
- The problems of having a too short second period is not as serious as those of having a too short first period, because after cutting off the presamples, we might not have enough samples left to draw a reliable result.
