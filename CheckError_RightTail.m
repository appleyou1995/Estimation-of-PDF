function y = CheckError_RightTail(mu, sigma, k, ...
                                  K_R0, K_R1, ...
                                  EMP_CDF_R0, EMP_PDF_R0, EMP_PDF_R1)
y = nan(1, 3);                                                             % Construct Space

% CDF at R0
z = (K_R0 - mu) / sigma;
y(1) = - (1 + k * z)^(- 1 / k) - log(EMP_CDF_R0);       
clear z 

% PDF at R0
z = (K_R0 - mu) / sigma;
y(2) = log((1 + k * z)^(- 1 - 1 / k)) - (1 + k * z)^(- 1 / k) - log(sigma) - log(EMP_PDF_R0);  
clear z 

% PDF at R1
z = (K_R1 - mu) / sigma;
y(3) = log((1 + k * z)^(- 1 - 1 / k)) - (1 + k * z)^(- 1 / k) - log(sigma) - log(EMP_PDF_R1); 
clear z    
end                                       