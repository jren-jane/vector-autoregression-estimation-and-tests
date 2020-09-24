%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% relative rejection frequencies %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% (01/896410) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ----------------------------------------------------------------------- %
%                                                                         %
%                           Introduction                                  %
%          --------------------------------------------------             %
%                                                                         %
%   This is the main file for computing the relative rejection            %
%   frequencies for the BP and SS test based on p-values obtained from    %
%   the asymptotic chi square distribution as well as on the bootstrap    %
%   p-values. The main file is based on a function for chowtest, a        %
%   function for Monte Carlo simulation, and a function for AIC           %
%   information criteria.                                                 %
%                                                                         %
% ----------------------------------------------------------------------- %

% ----------------------------------------------------------------------- %
%                         Declaration of Variables                        %
%          --------------------------------------------------             %
%                                                                         %
% (1)   Pmax: the maximum number of lags                                  %
%                                                                         %
% (2)   M: the number of Monte Carlo simulations                          %
%                                                                         %
% (3)   s: the number of periods to be discarded                          %
%                                                                         %
% (4)   indic: indicator of whether to include an intercept               %                                     
%                                                                         %
% (5)   siglvl: significance level                                        %
%                                                                         %
% (6)   T: length of data                                                 %
%                                                                         %
% (7)   Tb: break point                                                   %
%                                                                         %
% (8)   Yrep: bootstrapped time series returned by the other function     %
%       "MonteCarlo.m"                                                    %
%                                                                         %
% (9)   test: a matrix which temporarily retrieves the data generated in  %
%       the h-th replication                                              %
%                                                                         %
% (10)  pval_MC_bp_chisquare: a column vector storing p-values from the   %
%       BP test based on asymptotic chi-square distribution               % 
%                                                                         %
% (11)  pval_MC_ss_chisquare: a column vector storing p-values from the   %
%       SS test based on asymptotic chi-square distribution               %
%                                                                         %
% (12)  pval_MC_bp_bootstrap: a column vector storing p-values from the   %
%       BP test based on bootstrapped distribution                        %    
%                                                                         %
% (13)  pval_MC_ss_bootstrap: a column vector storing p-values from the   %
%       BP test based on bootstrapped distribution                        %
%                                                                         %
% (14)  rejfq_bp_chisquare_T1_Tb1: relative rejection frequency for the   %
%       BP test based on p-values obtained from asymptotic chi-square     %
%       distribution. The number of Monte Carlo simulation is 80 and the  %
%       break point is 0.5T.                                              %
%                                                                         %
% (15)  rejfq_ss_chisquare_T1_Tb1: relative rejection frequency for the   %
%       SS test based on p-values obtained from asymptotic chi-square     %
%       distribution. The number of Monte Carlo simulation is 80 and the  %
%       break point is 0.5T.                                              %
%                                                                         %
% (16)  rejfq_bp_bootstrap_T1_Tb1: relative rejection frequency for the   %
%       BP test based on p-values obtained from bootstrapped              %
%       distribution. The number of Monte Carlo simulation is 80 and the  %
%       break point is 0.5T.                                              %
%                                                                         %
% (17)  rejfq_ss_bootstrap_T1_Tb1: relative rejection frequency for the   %
%       SS test based on p-values obtained from bootstrapped              %
%       distribution. The number of Monte Carlo simulation is 80 and the  %
%       break point is 0.5T.                                              %
%                                                                         %
% (18)  rejfq_bp_chisquare_T2_Tb1: relative rejection frequency for the   %
%       BP test based on p-values obtained from asymptotic chi-square     %
%       distribution. The number of Monte Carlo simulation is 500 and the %
%       break point is 0.5T.                                              %
%                                                                         %
% (19)  rejfq_ss_chisquare_T2_Tb1: relative rejection frequency for the   %
%       SS test based on p-values obtained from asymptotic chi-square     %
%       distribution. The number of Monte Carlo simulation is 500 and the %
%       break point is 0.5T.                                              %
%                                                                         %
% (19)  rejfq_bp_bootstrap_T2_Tb1: relative rejection frequency for the   %
%       BP test based on p-values obtained from bootstrapped              %
%       distribution. The number of Monte Carlo simulation is 500 and the %
%       break point is 0.5T.                                              %
%                                                                         %
% (20)  rejfq_ss_bootstrap_T2_Tb1: relative rejection frequency for the   %
%       SS test based on p-values obtained from bootstrapped              %
%       distribution. The number of Monte Carlo simulation is 500 and the %
%       break point is 0.5T.                                              %
%                                                                         %
% (21)  rejfq_bp_chisquare_T1_Tb2: relative rejection frequency for the   %
%       BP test based on p-values obtained from asymptotic chi-square     %
%       distribution. The number of Monte Carlo simulation is 80 and the  %
%       break point is 0.2T.                                              %
%                                                                         %
% (22)  rejfq_ss_chisquare_T1_Tb2: relative rejection frequency for the   %
%       SS test based on p-values obtained from asymptotic chi-square     %
%       distribution. The number of Monte Carlo simulation is 80 and the  %
%       break point is 0.2T.                                              %
%                                                                         %
% (23)  rejfq_bp_bootstrap_T1_Tb2: relative rejection frequency for the   %
%       BP test based on p-values obtained from bootstrapped              %
%       distribution. The number of Monte Carlo simulation is 80 and the  %
%       break point is 0.2T.                                              %
%                                                                         %
% (24)  rejfq_ss_bootstrap_T1_Tb2: relative rejection frequency for the   %
%       SS test based on p-values obtained from bootstrapped              %
%       distribution. The number of Monte Carlo simulation is 80 and the  %
%       break point is 0.2T.                                              %
%                                                                         %
% (25)  rejfq_bp_chisquare_T2_Tb2: relative rejection frequency for the   %
%       BP test based on p-values obtained from asymptotic chi-square     %
%       distribution. The number of Monte Carlo simulation is 500 and the %
%       break point is 0.2T.                                              %
%                                                                         %
% (26)  rejfq_ss_chisquare_T2_Tb2: relative rejection frequency for the   %
%       SS test based on p-values obtained from asymptotic chi-square     %
%       distribution. The number of Monte Carlo simulation is 500 and the %
%       break point is 0.2T.                                              %
%                                                                         %
% (27)  rejfq_bp_bootstrap_T2_Tb2: relative rejection frequency for the   %
%       BP test based on p-values obtained from bootstrapped              %
%       distribution. The number of Monte Carlo simulation is 500 and the %
%       break point is 0.2T.                                              %
%                                                                         %
% (28)  rejfq_ss_bootstrap_T2_Tb2: relative rejection frequency for the   %
%       SS test based on p-values obtained from bootstrapped              %
%       distribution. The number of Monte Carlo simulation is 500 and the %
%       break point is 0.2T.                                              %
%                                                                         %
%                                                                         %
% ----------------------------------------------------------------------- %
 
% ----------------- Specify the number of simulations ------------------- %

    Pmax = 4;
    M = 500;
    s = 50;
    indic = 1;                                                              % The model has an intercept
    siglvl = 0.05;                                                          % the significance level  

% ---------------- experiment with T = 80 & Tb = 0.5T ------------------- %

    T = 80;
    Tb = 0.5 * T;

    pval_MC_bp_chisquare = zeros(M , 1);                                    % Initialize the vector for storing the p-values
    pval_MC_ss_chisquare = zeros(M , 1);
    pval_MC_bp_bootstrap = zeros(M , 1);
    pval_MC_ss_bootstrap = zeros(M , 1);

for h = 1 : M
    
    [Yrep] = MonteCarlo(M, s, T);
    
    test = Yrep(2 * h - 1 : 2 * h , :);
    
    p = find_lag_AIC(test , Pmax , indic);
    
    [pval_bp_chisquare , pval_ss_chisquare , pval_bp_bootstrap , pval_ss_bootstrap] = chowtest_chisquare_bootstrap(test , Tb , p , indic);
    pval_MC_bp_chisquare(h) = pval_bp_chisquare;
    pval_MC_ss_chisquare(h) = pval_ss_chisquare;
    pval_MC_bp_bootstrap(h) = pval_bp_bootstrap;
    pval_MC_ss_bootstrap(h) = pval_ss_bootstrap;
    
end

    rejfq_bp_chisquare_T1_Tb1 = sum(pval_MC_bp_chisquare < siglvl) / M      
    rejfq_ss_chisquare_T1_Tb1 = sum(pval_MC_ss_chisquare < siglvl) / M
    rejfq_bp_bootstrap_T1_Tb1 = sum(pval_MC_bp_bootstrap < siglvl) / M
    rejfq_ss_bootstrap_T1_Tb1 = sum(pval_MC_ss_bootstrap < siglvl) / M

% ------------------------ experiment with T = 500 & Tb = 0.5T ---------------------------- %

    T = 500;
    Tb = 0.5 * T;


for h = 1 : M
    
    [Yrep] = MonteCarlo(M, s, T);
    
    test = Yrep(2 * h - 1 : 2 * h , :);
    
    p = find_lag_AIC(test , Pmax , indic);
    
    [pval_bp_chisquare , pval_ss_chisquare , pval_bp_bootstrap , pval_ss_bootstrap] = chowtest_chisquare_bootstrap(test , Tb , p , indic);
    pval_MC_bp_chisquare(h) = pval_bp_chisquare;
    pval_MC_ss_chisquare(h) = pval_ss_chisquare;
    pval_MC_bp_bootstrap(h) = pval_bp_bootstrap;
    pval_MC_ss_bootstrap(h) = pval_ss_bootstrap;
    
end

    rejfq_bp_chisquare_T2_Tb1 = sum(pval_MC_bp_chisquare < siglvl) / M
    rejfq_ss_chisquare_T2_Tb1 = sum(pval_MC_ss_chisquare < siglvl) / M
    rejfq_bp_bootstrap_T2_Tb1 = sum(pval_MC_bp_bootstrap < siglvl) / M
    rejfq_ss_bootstrap_T2_Tb1 = sum(pval_MC_ss_bootstrap < siglvl) / M

% ------------------------ experiment with T = 80 & Tb = 0.2T ---------------------------- %

    T = 80;
    Tb = 0.2 * T;

for h = 1 : M
    
    [Yrep] = MonteCarlo(M, s, T);
    
    test = Yrep(2 * h - 1 : 2 * h , :);
    
    p = find_lag_AIC(test , Pmax , indic);
    
    [pval_bp_chisquare , pval_ss_chisquare , pval_bp_bootstrap , pval_ss_bootstrap] = chowtest_chisquare_bootstrap(test , Tb , p , indic);
    pval_MC_bp_chisquare(h) = pval_bp_chisquare;
    pval_MC_ss_chisquare(h) = pval_ss_chisquare;
    pval_MC_bp_bootstrap(h) = pval_bp_bootstrap;
    pval_MC_ss_bootstrap(h) = pval_ss_bootstrap;
    
end

    rejfq_bp_chisquare_T1_Tb2 = sum(pval_MC_bp_chisquare < siglvl) / M  
    rejfq_ss_chisquare_T1_Tb2 = sum(pval_MC_ss_chisquare < siglvl) / M
    rejfq_bp_bootstrap_T1_Tb2 = sum(pval_MC_bp_bootstrap < siglvl) / M
    rejfq_ss_bootstrap_T1_Tb2 = sum(pval_MC_ss_bootstrap < siglvl) / M

% ------------------------ experiment with T = 500 & Tb = 0.2T ---------------------------- %

    T = 500;
    Tb = 0.2 * T;

for h = 1 : M
    
    [Yrep] = MonteCarlo(M, s, T);
    
    test = Yrep(2 * h - 1 : 2 * h , :);
    
    p = find_lag_AIC(test , Pmax , indic);
    
    [pval_bp_chisquare , pval_ss_chisquare , pval_bp_bootstrap , pval_ss_bootstrap] = chowtest_chisquare_bootstrap(test , Tb , p , indic);
    pval_MC_bp_chisquare(h) = pval_bp_chisquare;
    pval_MC_ss_chisquare(h) = pval_ss_chisquare;
    pval_MC_bp_bootstrap(h) = pval_bp_bootstrap;
    pval_MC_ss_bootstrap(h) = pval_ss_bootstrap;
    
end

    rejfq_bp_chisquare_T2_Tb2 = sum(pval_MC_bp_chisquare < siglvl) / M
    rejfq_ss_chisquare_T2_Tb2 = sum(pval_MC_ss_chisquare < siglvl) / M
    rejfq_bp_bootstrap_T2_Tb2 = sum(pval_MC_bp_bootstrap < siglvl) / M
    rejfq_ss_bootstrap_T2_Tb2 = sum(pval_MC_ss_bootstrap < siglvl) / M