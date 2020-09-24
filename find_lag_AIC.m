%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Function for AIC information criterion  %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%   (01/896410)   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ----------------------------------------------------------------------- %
%                                                                         %
%                           Introduction                                  %
%         --------------------------------------------------              %
%                                                                         %
%   This function determines which order should be chosen based on AIC    %
%   information criterion. The inputs of the function are the data, the    %
%   maximum order up until which we compare, and the indicator which      %
%   implies whether the intercept has been included. The output of the    %
%   function is the optimal choice of order by AIC.                       %
%                                                                         %
% ----------------------------------------------------------------------- %

% ----------------------------------------------------------------------- %
%                         Declaration of Variables                        %
%          --------------------------------------------------             %
%                                                                         %
% (1)   Data: the data input, in our problem a given (T+p)*K matrix       %
%                                                                         %
% (2)   Pmax: the maximum number of lags                                  %
%                                                                         %
% (3)   indic: indicator of whether to include an intercept               %
%                                                                         %
% (4)   AICest: the lag order suggested by the AIC information criterion  %
%                                                                         %
% (5)   T: length of data                                                 %
%                                                                         %
% (6)   K: the dimention of a single observation                          %
%                                                                         %
% (7)   AIC: a column vector which stores the AIC criterion of all lag    %
%       orders                                                            %
%                                                                         %
% ----------------------------------------------------------------------- %


    function [AICest] = find_lag_AIC(Data , Pmax , indic)

% --------------------------- Initial Set-up ---------------------------- %

    T      =   size(Data , 2) - Pmax;                                          % The number of observations except pre-sample.
    K      =   size(Data , 1);                                                 % The dimension of the single observation.

    Y      =   Data(:, (Pmax + 1) : end);                                      % Collecting the remaining observations to matrix Y.

    AIC    =   zeros(Pmax + 1 , 1);

% -------------- Arrange the Data to the Compact Form ------------------- %

    for m = 0 : Pmax

    Z  = zeros(K * m , T);

    if indic  ==  1;
        for i  =  0 : (m - 1)
    Z(1 + i * K : (i + 1) * K , :) = Data(: , Pmax - i : T + Pmax - 1 - i);
        end
    Z  = [ones(1 , T) ; Z];
    else
        for i  =  0 : (m - 1)
    Z(1 + i * K : (i + 1) * K , :) = Data(: , Pmax - i : T + Pmax - 1 - i);
        end
    end

% -------------- LS Estimator of the VAR(p) Model  --------------------- %

    BE     =   (Y * Z') / (Z * Z');

% ---------- Calculating Residual Matrix in Compact Form --------------- %

    UE     =   Y - BE * Z;

% ---- Calculating the unadjusted covariance matrix of the residuals --- %
                                  
    Utilde     =   1 / T * (UE * UE');
                                                          
% ---------------------------AIC criterion------------------------------ %

    AIChelp  =  log(det(Utilde)) + 2 * m * K^2 / T;
    AIC(m + 1) = AIChelp;
   
    end
                        
      [ ~ , lagest]   =   min(AIC);
      AICest = lagest-1;
end
