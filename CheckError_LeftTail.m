function y = CheckError_LeftTail(mu, sigma, k, ...
                                 K_L0, K_L1, ...
                                 EMP_CDF_L0, EMP_PDF_L0, EMP_PDF_L1)
y = nan(1, 3);                                                             % Construct Space

% CDF at L0
z = (K_L0 - mu) / sigma;
y(1) = - (1 + k * z)^(- 1 / k) - log(1 - EMP_CDF_L0);       
clear z 

% PDF at L0
z = (K_L0 - mu) / sigma;
y(2) = log((1 + k * z)^(- 1 - 1 / k)) - (1 + k * z)^(- 1 / k) - log(sigma) - log(EMP_PDF_L0);  
clear z 

% PDF at L1
z = (K_L1 - mu) / sigma;
y(3) = log((1 + k * z)^(- 1 - 1 / k)) - (1 + k * z)^(- 1 / k) - log(sigma) - log(EMP_PDF_L1); 
clear z    
end          