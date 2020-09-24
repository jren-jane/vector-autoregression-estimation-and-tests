%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Function for Chow Test  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%       (01/896410)       %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ----------------------------------------------------------------------- %
%                                                                         %
%                           Introduction                                  %
%         --------------------------------------------------              %
%                                                                         %
%   This function determines which order should be chosen based on AIC    %
%   and SC information criteria, respectively. The inputs of the function %
%   are the data, the maximum order up until which we compare, and the    %
%   indicator which implies whether the intercept has been included. The  %
%   output of the function is the optimal choice of order by AIC and SC,  %
%   respectively.                                                         %
%                                                                         %
% ----------------------------------------------------------------------- %

% ----------------------------------------------------------------------- %
%                         Declaration of Variables                        %
%          --------------------------------------------------             %
%                                                                         %
% (1）Data: the data input, in our problem a given (T+p)*K matrix         %
%                                                                         %
% (2) Tb: the break point                                                 %
%                                                                         %
% (3) p: number of lags                                                   %
%                                                                         %
% (4) indic: indicator of whether to include an intercept                 %
%                                                                         %
%（5）pval_bp_chisquare: the p-value of the breakpoint test based on an    %
%    asymptotic chisquare distribution                                    %
%                                                                         %
% (6)pval_ss_chisquare: the p-value of the sample split test based on an  %
%    asymptotic chisquare distribution                                    %
%                                                                         %
% (7)pval_bp_bootstrap:the p-value of the breakpoint test based on a      %
%    bootstrapped distribution                                            %
%                                                                         %
% (8)pval_ss_bootstrap:the p-value of the sample split test based on a    %
%    bootstrapped distribution                                            %
% ----------------------------------------------------------------------- %

function [pval_bp_chisquare , pval_ss_chisquare , pval_bp_bootstrap , pval_ss_bootstrap] = chowtest_chisquare_bootstrap(Data , Tb , p , indic)

% -------- test based on the asymptotic chi-square distribution --------- %

% --------------------------- Initial Set-up ---------------------------- %

    K      =    size(Data , 1);                                             % the dimension of the single observation
    T      =    size(Data , 2) - p;                                         % the number of all observations except pre-samples
    T1     =    Tb - 1;                                                     % the number of observations in period 1 except pre-sample
    T2     =    T - T1 - p;                                                 % the number of observations in period 2 except pre-sample

    Y1     =    Data(: , (p + 1) : (p + T1));                               % observations in period 1 except pre-samples.
    Y2     =    Data(: , (p + T1 + p + 1) : (p + T1 + p + T2));             % observations in period 2 except pre-samples.
    Psam   =    Data(: , 1 : p);                                            % pre-samples are extracted from the data.
    Y      =    Data(: , (p + 1) : end);                                    % all observations except pre-sample

% -------------- Arrange the Data to the Compact Form ------------------- %

    Z1  =   zeros(K * p , T1);                                              % compact form before the break
      if indic == 1
        for i = 0 : (p - 1)
          Z1(1 + i * K : (i + 1) * K, :) = Data(: , p - i : T1 + (p - i) - 1);
        end
      Z1  =  [ones(1 , T1) ; Z1];
      else
        for i = 0 : (p - 1)
          Z1(1 + i * K : (i + 1) * K, :) = Data(: , p - i : T1 + (p - i) - 1);
        end
      end

    Z2  =   zeros(K * p , T2);                                              % compact form after the break
      if indic == 1;
        for i = 0 : (p - 1)
          Z2(1 + i * K : (i + 1) * K, :) = Data(: , (T1 + p) + (p - i) : T2 + (T1 + p) + (p - i) - 1);
        end
      Z2  =   [ones(1 , T2) ; Z2];
      else
        for i = 0 : (p - 1)
          Z2(1 + i * K : (i + 1) * K, :) = Data(: , (Tb + p - 1) + (p - i) : T2 + (Tb + p - 1) + (p - i) - 1);
        end
      end

    Z   =   zeros(K * p , T);                                               % compact form for the full model
      if indic  ==  1
        for i  =  0 : (p - 1)
          Z(1 + i * K : (i + 1) * K , :) = Data(: , p - i : T + p - 1 - i);
        end
      Z   =  [ones(1 , T) ; Z];
      else
        for i  =  0 : (p - 1)
          Z(1 + i * K : (i + 1) * K , :) = Data(: , p - i : T + p - 1 - i);
        end
      end

% -------------- LS Estimator of the VAR(p) Model  --------------------- %

    BE1    =   (Y1 * Z1') / (Z1 * Z1');
    BE2    =   (Y2 * Z2') / (Z2 * Z2');
    BE     =   (Y * Z') / (Z * Z');

% ------------------ Residual Matrix in Compact Form ------------------- %

    UE1    =   Y1 - BE1 * Z1;
    UE2    =   Y2 - BE2 * Z2;
    UE     =   Y - BE * Z;

% ----------------- Covariance Matrix of the Residuals ----------------- %

    sigma_tilde1   =   1 /  T1 * (UE1 * UE1');
    sigma_tilde2   =   1 /  T2 * (UE2 * UE2');
    sigma_tilde12  =   1 / (T1 + T2) * (UE(: , 1 : T1) * UE(: , 1 : T1)'+ UE(: , (T - T2 + 1) : end) * UE(: , (T - T2 + 1) : end)');
                                                            
% ---------------------- the Chow Test Statistics ---------------------- %
                                                            
    lambda_bp  =  (T1 + T2) * log( det( sigma_tilde12 ) ) - T1 * log( det(sigma_tilde1) ) - T2 * log ( det(sigma_tilde2) );  % the break point test statistic
    lambda_ss  =  (T1 + T2) * (log( det( sigma_tilde12 ) ) - log( det( (T1 * sigma_tilde1 + T2 * sigma_tilde2)/(T1 + T2)))); % the sample split test statistic
    
    lambda_bp_original=lambda_bp;                                           % Store the original break point test statistic to be later compared to the bootstrapped values
    lambda_ss_original=lambda_ss;                                           % Store the original sample split test statistic to be later compared to the bootstrapped values

% --------------------------- the p-values ----------------------------- %
                              
    dof_bp = K^2 * p + K + K * (K + 1) / 2;
    dof_ss = K^2 * p + K;
                                                            
    pval_bp_chisquare = 1 - chi2cdf(lambda_bp , dof_bp);
    pval_ss_chisquare = 1 - chi2cdf(lambda_ss , dof_ss);
                              
  
    
    
% -------- test based on the distribution obtained by bootstraps--------- %

% --------------------- the number of bootstraps ------------------------ %

   B = 199;

% ----------------- Center the Residuals around their Means ------------- %
                                             
   UE_C   =   UE - repmat(mean(UE , 2) , 1 , T);
   
% --------- Draw Pseudorandom Integers based on Uniform Distr. ---------- %

   INDEX  =  randi([1 , T] , B , T);                                        % a B by T matrix which restores random integers
                                                                            % following a uniform distribution on interval [1 , T]
   
% ------------------- Initialize the test statistics -------------------- %
   
   lambda_bp_bootstrap = zeros(B , 1);                                      % This column vector stores the test stats from all replications
   lambda_ss_bootstrap = zeros(B , 1);                                      % This column vector stores the test stats from all replications

% -----------------------------Bootstrapping----------------------------- %

     for  j    =  1 : B                                                     % j denotes the number of replication
       Data2    =    [Psam  zeros(K , T)];                                  % This matrix temporarily stores the data for current replication.
                                                                            % The initial values are obtained from the pre-sample values.

     for  i  =  1 : T                                                       % i denotes the number of period, excluding the pre-sample
         y_temp     =    rot90(Data2(: , (i : (i + p - 1))))';              % y_temp denotes ystar(i-1)...ystar(i-p) listed in a column.
         y_temp     =    y_temp(:);                                         % In these two lines, a K*p matrix is transformed into a column vector
       if indic == 1
         y_temp   =    [1 ; y_temp];                                        % If there is an intercept, 1 should be added to the top of the column y_temp
       end
     Data2(: , i + p) = BE * y_temp + UE_C(: , INDEX(j , i));               % Fill in ystar(i) based on ystar(i-1)...ystar(i-p).
                                                                            % INDEX(j , i) denotes the index for the draw in the vector UE_C
                                                                            % in the ith period and the jth replication
     end

% -------------------------- Initial Set-up ----------------------------- %
                                                                            % The same estimation process as before is applied to each bootstrap.
    K      =    size(Data2 , 1);                                            % The dimension of the single observation.
    T      =    size(Data2 , 2) - p;                                        % The number of all observations except pre-samples.
    T1     =    Tb - 1;                                                     % The number of observations in period 1 except pre-sample.
    T2     =    T - T1 - p;                                                 % The number of observations in period 2 except pre-sample.
                              
    Y1     =    Data2(: , (p + 1) : (p + T1));                              % observations in period 1 except pre-samples.
    Y2     =    Data2(: , (p + T1 + p + 1) : (p + T1 + p + T2));            % observations in period 2 except pre-samples.
    Y      =    Data2(: , (p + 1) : end);                                   % all observations except pre-sample
                              
% -------------- Arrange the Data to the Compact Form ------------------- %
                              
    Z1  =   zeros(K * p , T1);
      if indic == 1;
        for i = 0 : (p - 1)
          Z1(1 + i * K : (i + 1) * K, :) = Data(: , p - i : T1 + (p - i) - 1);
        end
      Z1  =  [ones(1 , T1) ; Z1];
      else
        for i = 0 : (p - 1)
          Z1(1 + i * K : (i + 1) * K, :) = Data(: , p - i : T1 + (p - i) - 1);
        end
      end
                              
    Z2  =   zeros(K * p , T2);
      if indic == 1;
        for i = 0 : (p - 1)
          Z2(1 + i * K : (i + 1) * K, :) = Data(: , (T1 + p) + (p - i) : T2 + (T1 + p) + (p - i) - 1);
        end
      Z2  =   [ones(1 , T2) ; Z2];
      else
        for i = 0 : (p - 1)
          Z2(1 + i * K : (i + 1) * K, :) = Data(: , (Tb + p - 1) + (p - i) : T2 + (Tb + p - 1) + (p - i) - 1);
        end
      end
                              
    Z   =   zeros(K * p , T);
      if indic  ==  1;
        for i  =  0 : (p - 1)
          Z(1 + i * K : (i + 1) * K , :) = Data(: , p - i : T + p - 1 - i);
        end
      Z   =  [ones(1 , T) ; Z];
      else
        for i  =  0 : (p - 1)
          Z(1 + i * K : (i + 1) * K , :) = Data(: , p - i : T + p - 1 - i);
        end
      end
                              
% ---------------- LS Estimator of the VAR(p) Model  ------------------- %
                              
    BE1    =   (Y1 * Z1') / (Z1 * Z1');
    BE2    =   (Y2 * Z2') / (Z2 * Z2');
    BE     =   (Y * Z') / (Z * Z');
                              
% ------------------ Residual Matrix in Compact Form ------------------- %
                              
    UE1    =   Y1 - BE1 * Z1;
    UE2    =   Y2 - BE2 * Z2;
    UE     =   Y - BE * Z;
                              
% ----------------- Covariance Matrix of the Residuals ----------------- %
                              
    sigma_tilde1   =   1 /  T1 * (UE1 * UE1');
    sigma_tilde2   =   1 /  T2 * (UE2 * UE2');
    sigma_tilde12  =   1 / (T1 + T2) * (UE(: , 1 : T1) * UE(: , 1 : T1)' + UE(: , (T - T2 + 1) : end) * UE(: , (T - T2 + 1) : end)');
                                                                                          
% ----------------------- the Chow Test Statistic ---------------------- %
                                                                                          
    lambda_bp  =  (T1 + T2) * log( det( sigma_tilde12 ) ) - T1 * log( det(sigma_tilde1) ) - T2 * log ( det(sigma_tilde2) );
    lambda_ss  =  (T1 + T2) * (log( det( sigma_tilde12 ) ) - log( det( (T1 * sigma_tilde1 + T2 * sigma_tilde2)/(T1 + T2))));
                              
                              
% ----------Fill in the series of test stats from replications---------- %
     
    lambda_bp_bootstrap(j , :) = lambda_bp;                                 % Fill in the break point test stat from the current replication
    lambda_ss_bootstrap(j , :) = lambda_ss;                                 % Fill in the sample split test stat from the current replication
     
     end                                                                    % End the for loop after collecting all the test stats from bootstraps.

% ------------------------Compute the p-values-------------------------- %

   pval_bp_bootstrap = sum(lambda_bp_bootstrap >= lambda_bp_original)/B;    % Calculate the number of bootstrapped test stats which are larger than the original test stat
                                                                            % and divide it by the total # of replications to obtain the p-value for the break point test

   pval_ss_bootstrap = sum(lambda_ss_bootstrap >= lambda_ss_original)/B;    % Calculate the number of bootstrapped test stats which are larger than the original test stat
                                                                            % and divide it by the total # of replications to obtain the p-value for the sample split test

   end
