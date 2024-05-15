clear; clc


%% Setting

% Specific Date
Target_Date = 19960717;

% Specific Time-to-Maturity (LB)
Target_TTM = 29;


%% Load Data

% [1.  SecID              | 2.  Date (YYYYMMDD)    | 3.  TTM (Days)    | 4.  CPflag | 5.  K  | 6. S | 7. F | 
% [8.  Option Price (Bid) | 9.  Option Price (Ask) | 10. Open Interest | 11. Volume | 12. IV | 
% [13. Delta              | 14. Gamma              | 15. Theta         | 16. Vega                          ]

FileName = ['OP' num2str(fix(Target_Date / 10000)) '_' num2str(fix(rem(Target_Date, 10000) / 100)) '.txt'];
Data = load(['D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\CDI Method\Data\99 姿穎學姊提供\20240417\IndexOptions19962019_SP500\' FileName]);
clear FileName


% Index of Data
Index_ID = 1;
Index_Date = 2;
Index_TTM = 3;
Index_CPFlag = 4;
Index_K = 5;
Index_S = 6;
Index_F = 7;                                                               % Forward price (Theoretical)
Index_OP_Bid = 8;
Index_OP_Ask = 9;
Index_OI = 10;
Index_V = 11;
Index_IV = 12;
Index_Delta = 13;
Index_Gamma = 14;
Index_Theta = 15;
Index_Vega = 16;
Index_RF = Index_Vega + 1;                                                 % Construction
Index_DY = Index_Vega + 2; 


% Specific Data
Index = find((Data(:, Index_Date)==Target_Date) & ...
             (Data(:, Index_TTM) >= Target_TTM));
Data = Data(Index, :);                                                     % Update: Data
clear Index

% Specific Time-to-Maturity 
Target_TTM = min(Data(:, Index_TTM));                                      % Update: Target_TTM

Index = find(Data(:, Index_TTM)==Target_TTM);
Data = Data(Index, :);                                                     % Update: Data
clear Index


%% Correction of Expiration Date (Time-to-Maturity)

% Saturday to Friday
Date_EXP_WeekDay = weekday(datenum(num2str(Data(:, Index_Date)), 'yyyymmdd') ...
                           + Data(:, Index_TTM));      
Data(Date_EXP_WeekDay==7, Index_TTM) = Data(Date_EXP_WeekDay==7, Index_TTM) - 1;   % Update: Data    
clear Date_EXP_WeekDay            
    
% AM Settlement
% Calculation of Risk-Free Rate and Implied Volatility (IvyDB Reference Manual) 
Data(:, Index_TTM) = Data(:, Index_TTM) - 1;                               % Update: Data   

Target_TTM = Data(1, Index_TTM);


%% Data Filtering

% Bid >= 0.375
Index = find(Data(:, Index_OP_Bid) >= 0.375);
Data = Data(Index, :);                                                     % Update: Data
clear Index


% Out-of-the-Money (OTM) and Around At-the-Money (ATM) Options
Index = find((Data(:, Index_CPFlag)==1) & ...
             (Data(:, Index_K) >= (Data(:, Index_S))) & ...
             (Data(:, Index_K) <= (Data(:, Index_S))));
Data(Index, Index_CPFlag) = 31;                                            % Update: Data (CP Flag)           
clear Index

Index = find((Data(:, Index_CPFlag)==2) & ...
             (Data(:, Index_K) >= (Data(:, Index_S))) & ...
             (Data(:, Index_K) <= (Data(:, Index_S))));
Data(Index, Index_CPFlag) = 32;                                            % Update: Data (CP Flag)           
clear Index
          
Index = find(((Data(:, Index_CPFlag)==1) & (Data(:, Index_K) >= Data(:, Index_S))) | ...
             ((Data(:, Index_CPFlag)==2) & (Data(:, Index_K) <= Data(:, Index_S))) | ...
             (fix(Data(:, Index_CPFlag) / 10)==3));      
Data = Data(Index, :);                                                     % Update: Data 
clear Index 


%% Calculation of Dividend Yield

% [1. SecID | 2. Date (YYYYMMDD) | 3. Dividend Yield (Annualized)]

FileName = ['IndexDivYield19962019.txt'];
Data_DY = load(['D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\CDI Method\Data\99 姿穎學姊提供\20240417\' FileName]);
clear FileName

Data(:, Index_DY) = Data_DY(Data_DY(:, Index_Date)==Target_Date, end);     % Annualized
% clear Data_DY


%% Calculation of Risk-Free Rate

% [1. Date (YYYYMMDD) | 2. TTM (Days) | 3. Risk-Free Rate (Annualized)]

FileName = ['RiskFreeRate19962019.txt'];
Data_RF = load(['D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\CDI Method\Data\99 姿穎學姊提供\20240417\' FileName]);
clear FileName

% RF_Target = RF_TTM(Data, Date_Target, TTM_Target)
Data(:, Index_RF) = RF_TTM(Data_RF, Data(:, Index_Date), Data(:, Index_TTM));   % Annualized
% clear Data_RF


%% Determine the Range of Density Function

S0 = Data(1, Index_S);
TTM = Data(1, Index_TTM) / 365;                                            % Annualized
DY = Data(1, Index_DY);                                                    % Annualized
RF = Data(1, Index_RF);                                                    % Annualized

% All Possible Range
K_Low = (0.3 / 100) * S0;                                                  % Hollstein and Prokopczuk (2016 JFQA) 
K_High = 3 * S0;                                                           % Hollstein and Prokopczuk (2016 JFQA) 
Smooth_AllK = K_Low:0.05:K_High;

% clear K_Low K_High


%% Spline Approximation (Empirical)

% Fourth-Order Spline Approximation  
LSS = spap2(1, 4, Data(:, Index_K), Data(:, Index_IV));


% Smoothed Implied Volatility
Smooth_K = Smooth_AllK((Smooth_AllK >= min(Data(:, Index_K))) & ...
                       (Smooth_AllK <= max(Data(:, Index_K))));
Smooth_IV = fnval(LSS, Smooth_K);
% clear LSS


% Smoothed Option Price
S0 = S0 * ones(size(Smooth_K));                                            % Update: S0
TTM = TTM * ones(size(Smooth_K));                                          % Update: TTM
DY = DY * ones(size(Smooth_K));                                            % Update: DY
RF = RF * ones(size(Smooth_K));                                            % Update: RF

Smooth_OP = blsprice(S0, Smooth_K, RF, TTM, Smooth_IV, DY); 

% clear Smooth_IV
% clear S0 DY


%% Plot Figure: Smoothed Implied Volatility (Fourth Order Spline Approximation)

figure

plot(Smooth_K, Smooth_IV, ...
     'LineWidth', 2, ...
     'Color', 'r');   
hold on

plot(Data(:, Index_K), Data(:, Index_IV), ...
     'LineStyle', 'none', ...     
     'LineWidth', 2, ...
     'Marker', 'o', ...     
     'MarkerSize', 6, ...
     'Color', 'b'); 
grid on

h1 = title(['Implied Volatilities with Interpolation and Smoothing (Date = ' datestr(datenum(num2str(Target_Date), 'yyyymmdd'), 'mmm dd, yyyy') ', Time to Maturity = ' num2str(Target_TTM) ' Days)']);
h2 = legend(['4th Degree Spline Interpolation (Without Knot)'], 'Raw Option Data');
h3 = xlabel('Strike Price');
h4 = ylabel('Implied Volatility');

set([h1 h2 h3 h4], ...
    'FontSize', 14, ...
    'FontWeight', 'bold', ...
    'LineWidth', 2)
clear h1 h2 h3 h4


%% Risk-Neutral Density and Distribution (Empirical)

% Smooth_EMP_PDF and Smooth_EMP_CDF (w.r.t. Smooth_K)
TTM = TTM(1);                                                              % Update: TTM
RF = RF(1); 

% Risk-Neutral Density (PDF)
% K_{2} to K_{N - 1}
Smooth_EMP_PDF = exp(RF * TTM) * ...
                 (Smooth_OP(3:end) - 2 * Smooth_OP(2:(end - 1)) + Smooth_OP(1:(end - 2))) ./ ...
                 ((Smooth_K(3:end) - Smooth_K(2:(end - 1))) .^ 2);

% Risk-Neutral Distribution (CDF)
% K_{2} to K_{N - 1}
Smooth_EMP_CDF = exp(RF * TTM) * ...
                 ((Smooth_OP(3:end) - Smooth_OP(1:(end - 2))) ./ ...
                 (Smooth_K(3:end) - Smooth_K(1:(end - 2)))) + 1;
% clear Smooth_OP
% clear RF TTM

% K_{2} to K_{N - 1}
Smooth_K = Smooth_K(2:(end - 1));                                          % Update: Smooth_K
Smooth_IV = Smooth_IV(2:(end - 1));                                        % Update: Smooth_IV (Just Record)
Smooth_OP = Smooth_OP(2:(end - 1));


%% Plot Figure: CDF

figure

plot(Smooth_K, Smooth_EMP_CDF, ...
     'LineWidth', 2, ...
     'Color', 'r');  
hold on

set(gca,'XTick', min(Smooth_AllK):10:max(Smooth_AllK));
xlim([540 690]);
grid on

h1 = title(['CDF (Date = ' datestr(datenum(num2str(Target_Date), 'yyyymmdd'), 'mmm dd, yyyy') ', Time to Maturity = ' num2str(Target_TTM) ' Days)']);
h2 = legend('Empirical CDF', 'Location', 'northwest');
h3 = xlabel('Level of S&P 500 Index');
h4 = ylabel('Cummulated Density');

set([h1 h2 h3 h4], ...
    'FontSize', 14, ...
    'FontWeight', 'bold', ...
    'LineWidth', 2)
clear h1 h2 h3 h4


%% Risk-Neutral Density (Right-Tail Connection)

% Setting
BP_R0 = 0.92; 
BP_R1 = 0.95; 
 
% (BP_R1, K_R1, EMP_PDF_R1)
Index = find((Smooth_EMP_CDF - BP_R1) >= 0);
if length(Index) > 0
    BP_R1 = Smooth_EMP_CDF(Index(1));                                      % Update: BP_R1    
    
    K_R1 = Smooth_K(Index(1)); 
    EMP_PDF_R1 = Smooth_EMP_PDF(Index(1));
else
    BP_R1 = Smooth_EMP_CDF(end);                                           % Update: BP_R1    
    
    K_R1 = Smooth_K(end);     
    EMP_PDF_R1 = Smooth_EMP_PDF(end);      
end
BP_R0 = BP_R1 - 0.03;                                                  % Update: BP_R0 
clear Index 


% (BP_R0, K_R0, EMP_PDF_R0, EMP_CDF_R0)
Index = find((Smooth_EMP_CDF - BP_R0) >= 0);
BP_R0 = Smooth_EMP_CDF(Index(1));                                          % Update: BP_R0
K_R0 = Smooth_K(Index(1));
EMP_PDF_R0 = Smooth_EMP_PDF(Index(1));
EMP_CDF_R0 = Smooth_EMP_CDF(Index(1));
clear Index 

% Solve System of Equations (Three Conditions)
% (mu, sigma, k) = (Location, Scale, Shape)
syms x mu sigma k
z = (x - mu) / sigma;
eq1 = subs(exp(- (1 + k * z)^(- 1 / k)), 'x', K_R0) - EMP_CDF_R0;                                         % CDF at R0
eq2 = subs(((1 + k * z)^(- 1 - 1 / k)) * exp(- (1 + k * z)^(- 1 / k)) / sigma, 'x', K_R0) - EMP_PDF_R0;   % PDF at R0
eq3 = subs(((1 + k * z)^(- 1 - 1 / k)) * exp(- (1 + k * z)^(- 1 / k)) / sigma, 'x', K_R1) - EMP_PDF_R1;   % PDF at R1   

sol = solve([eq1, eq2, eq3], mu, sigma, k);


if (length(sol.mu)==0) | (length(sol.sigma)==0) | (length(sol.k)==0)
    clear eq1 eq2 eq3 sol
    eq1 = subs(- (1 + k * z)^(- 1 / k), 'x', K_R0) - log(EMP_CDF_R0);                                               % CDF at R0
    eq2 = subs(log((1 + k * z)^(- 1 - 1 / k)) - (1 + k * z)^(- 1 / k) - log(sigma), 'x', K_R0) - log(EMP_PDF_R0);   % PDF at R0
    eq3 = subs(log((1 + k * z)^(- 1 - 1 / k)) - (1 + k * z)^(- 1 / k) - log(sigma), 'x', K_R1) - log(EMP_PDF_R1);   % PDF at R1    
    
    sol = solve([eq1, eq2, eq3], mu, sigma, k);
else
end
% clear EMP_CDF_R0 EMP_PDF_R0 EMP_PDF_R1
clear x z eq1 eq2 eq3 


% (Parameters_GEV_R, Smooth_GEV_R_PDF, Smooth_GEV_R_CDF)
try
    mu = real(double(sol.mu));                                             % Update: mu (System to Value)
    sigma = real(double(sol.sigma));                                       % Update: sigma (System to Value)
    k = real(double(sol.k));                                               % Update: k (System to Value)

    % Right Tail of Risk-Neutral Density (PDF)
    Smooth_GEV_R_PDF = gevpdf(Smooth_AllK, ...
                              k, sigma, mu);
                      
    % Right Tail of Risk-Neutral Distribution (CDF)
    Smooth_GEV_R_CDF = gevcdf(Smooth_AllK, ...
                              k, sigma, mu);     
                          
    % Parameters 
    Parameters_GEV_R = [mu sigma k];     
    
    % Error of Three Conditions
    FitError_GEV_R = CheckError_RightTail(mu, sigma, k, ...
                                          K_R0, K_R1, ...
                                          EMP_CDF_R0, EMP_PDF_R0, EMP_PDF_R1);          
catch
end     
clear EMP_CDF_R0 EMP_PDF_R0 EMP_PDF_R1
clear mu sigma k sol


%% Risk-Neutral Density (Left-Tail Connection, Reverse Left to Right)

% Setting
BP_L0 = 0.05;
BP_L1 = 0.02;

for n = 1:10000
    n
    
    % Setting
    if n > 1
        BP_L1 = min(Smooth_EMP_CDF(Smooth_EMP_CDF > BP_L1));               % Update: BP_L1  
        BP_L0 = BP_L1 + 0.03;                                              % Update: BP_L0 
    else
    end
        
    % (BP_L1, K_L1, EMP_PDF_L1)
    Index = find((BP_L1 - Smooth_EMP_CDF) >= 0);
    if length(Index) > 0
        BP_L1 = Smooth_EMP_CDF(Index(end));                                % Update: BP_L1 
    
        K_L1 = Smooth_K(Index(end));  
        EMP_PDF_L1 = Smooth_EMP_PDF(Index(end));
    else
        BP_L1 = Smooth_EMP_CDF(1); % Update: BP_L1    
    
        K_L1 = Smooth_K(1);  
        EMP_PDF_L1 = Smooth_EMP_PDF(1);    
    end
    BP_L0 = BP_L1 + 0.03;                                                  % Update: BP_L0
    clear Index 

    % (BP_L0, K_L0, EMP_PDF_L0, EMP_CDF_L0)
    Index = find((BP_L0 - Smooth_EMP_CDF) >= 0);
    BP_L0 = Smooth_EMP_CDF(Index(end));                                    % Update: BP_L0  
    K_L0 = Smooth_K(Index(end));
    EMP_PDF_L0 = Smooth_EMP_PDF(Index(end));
    EMP_CDF_L0 = Smooth_EMP_CDF(Index(end));
    clear Index 

    % Solve System of Equations (Three Conditions)
    % (mu, sigma, k) = (Location, Scale, Shape)
    syms x mu sigma k
    z = (x - mu) / sigma;
    eq1 = subs(- (1 + k * z)^(- 1 / k), 'x', - K_L0) - log(1 - EMP_CDF_L0);                                           % CDF at L0
    eq2 = subs(log((1 + k * z)^(- 1 - 1 / k)) - (1 + k * z)^(- 1 / k) - log(sigma), 'x', - K_L0) - log(EMP_PDF_L0);   % PDF at L0
    eq3 = subs(log((1 + k * z)^(- 1 - 1 / k)) - (1 + k * z)^(- 1 / k) - log(sigma), 'x', - K_L1) - log(EMP_PDF_L1);   % PDF at L1
  
    sol = solve([eq1, eq2, eq3], mu, sigma, k);

    if (length(sol.mu)==0) | (length(sol.sigma)==0) | (length(sol.k)==0)
        clear eq1 eq2 eq3 sol
        eq1 = subs(exp(- (1 + k * z)^(- 1 / k)), 'x', - K_L0) - (1 - EMP_CDF_L0);                                   % CDF at L0
        eq2 = subs(((1 + k * z)^(- 1 - 1 / k)) * exp(- (1 + k * z)^(- 1 / k)) / sigma, 'x', - K_L0) - EMP_PDF_L0;   % PDF at L0
        eq3 = subs(((1 + k * z)^(- 1 - 1 / k)) * exp(- (1 + k * z)^(- 1 / k)) / sigma, 'x', - K_L1) - EMP_PDF_L1;   % PDF at L1
    
        sol = solve([eq1, eq2, eq3], mu, sigma, k);
    else
    end
    % clear EMP_CDF_L0 EMP_PDF_L0 EMP_PDF_L1
    clear x z eq1 eq2 eq3 


    % (Parameters_GEV_L, Smooth_GEV_L_PDF, Smooth_GEV_L_CDF)
    try
        mu = real(double(sol.mu));                                         % Update: mu (System to Value)
        sigma = real(double(sol.sigma));                                   % Update: sigma (System to Value)
        k = real(double(sol.k));                                           % Update: k (System to Value)

        % Left Tail of Risk-Neutral Density (PDF)
        Smooth_GEV_L_PDF = gevpdf(- Smooth_AllK, ...
                                  k, sigma, mu);
                      
        % Left Tail of Risk-Neutral Distribution (CDF)
        Smooth_GEV_L_CDF = 1 - gevcdf(- Smooth_AllK, ...
                                      k, sigma, mu);   
                          
        % Parameters 
        Parameters_GEV_L = [(- mu) sigma k];      
    
        % Error of Three Conditions
        FitError_GEV_L = CheckError_LeftTail(mu, sigma, k, ...
                                             - K_L0, - K_L1, ...
                                             EMP_CDF_L0, EMP_PDF_L0, EMP_PDF_L1);       
                                             
        % Check Existence and Fitting Error     
        if sum(isnan(Smooth_GEV_L_PDF))==0
            if max(abs(FitError_GEV_L)) < 0.00001
                break
            else
            end
        else
        end                                                
    catch
    end    
end
clear EMP_CDF_L0 EMP_PDF_L0 EMP_PDF_L1  
clear mu sigma k sol


%% Plot Figure

figure

% Risk-Neutral Density (Empirical)
plot(Smooth_K, Smooth_EMP_PDF, ...
     'LineWidth', 2, ...
     'Color', 'r');     
hold on
 
% Generalized Extreme Value Function (GEV, Right Tail)
try
    plot(Smooth_AllK, Smooth_GEV_R_PDF, ...
         'LineWidth', 2, ...
         'LineStyle', '--', ...
         'Color', 'b');    
    hold on
catch
end  

% Generalized Extreme Value Function (GEV, Left Tail)
try
    plot(Smooth_AllK, Smooth_GEV_L_PDF, ...
         'LineWidth', 2, ...
         'LineStyle', ':', ...
         'Color', 'g');    
    hold on
catch
end

% Connection Points
plot(Smooth_K(Smooth_K==K_R0), Smooth_EMP_PDF(Smooth_K==K_R0), ...
     'LineWidth', 2, ...
     'Marker', 'o', ...
     'Color', 'r');    
text(Smooth_K(Smooth_K==K_R0), Smooth_EMP_PDF(Smooth_K==K_R0), ...
     [num2str(100 * roundn(BP_R0, - 4)) '% \rightarrow  '], ...
     'HorizontalAlignment', 'right', ...
     'VerticalAlignment', 'baseline', ...
     'FontSize', 14, ...
     'FontWeight', 'bold', ...
     'LineWidth', 2, ...
     'Color', 'r');        
hold on

plot(Smooth_K(Smooth_K==K_R1), Smooth_EMP_PDF(Smooth_K==K_R1), ...
     'LineWidth', 2, ...
     'Marker', 'o', ...
     'Color', 'r');    
text(Smooth_K(Smooth_K==K_R1), Smooth_EMP_PDF(Smooth_K==K_R1), ...
     [num2str(100 * roundn(BP_R1, - 4)) '% \rightarrow  '], ...
     'HorizontalAlignment', 'right', ...
     'VerticalAlignment', 'baseline', ...
     'FontSize', 14, ...
     'FontWeight', 'bold', ...
     'LineWidth', 2, ...
     'Color', 'r');       
hold on

plot(Smooth_K(Smooth_K==K_L0), Smooth_EMP_PDF(Smooth_K==K_L0), ...
     'LineWidth', 2, ...
     'Marker', 'o', ...
     'Color', 'r');    
text(Smooth_K(Smooth_K==K_L0), Smooth_EMP_PDF(Smooth_K==K_L0), ...
     [num2str(100 * roundn(BP_L0, - 4)) '% \rightarrow  '], ...
     'HorizontalAlignment', 'right', ...
     'VerticalAlignment', 'baseline', ...
     'FontSize', 14, ...
     'FontWeight', 'bold', ...
     'LineWidth', 2, ...
     'Color', 'r');     
hold on

plot(Smooth_K(Smooth_K==K_L1), Smooth_EMP_PDF(Smooth_K==K_L1), ...
     'LineWidth', 2, ...
     'Marker', 'o', ...
     'Color', 'r');  
text(Smooth_K(Smooth_K==K_L1), Smooth_EMP_PDF(Smooth_K==K_L1), ...
     [num2str(100 * roundn(BP_L1, - 4)) '% \rightarrow  '], ...
     'HorizontalAlignment', 'right', ...
     'VerticalAlignment', 'baseline', ...
     'FontSize', 14, ...
     'FontWeight', 'bold', ...
     'LineWidth', 2, ...
     'Color', 'r');   

set(gca,'XTick', min(Smooth_AllK):50:max(Smooth_AllK));
xlim([min(Smooth_AllK) max(Smooth_AllK)]);
grid on

h1 = title(['Risk-Neutral Density and Fitted GEV Tail Functions (Date = ' datestr(datenum(num2str(Target_Date), 'yyyymmdd'), 'mmm dd, yyyy') ', Time to Maturity = ' num2str(Target_TTM) ' Days)']);
h2 = legend('Empirical RND', 'Right Tail GEV Function', 'Left Tail GEV Function', 'Connection Points');
h3 = xlabel('Level of S&P 500 Index');
h4 = ylabel('Density');

set([h1 h2 h3 h4], ...
    'FontSize', 14, ...
    'FontWeight', 'bold', ...
    'LineWidth', 2)
clear h1 h2 h3 h4


%% Combination of Full Density 

Smooth_RND = nan(size(Smooth_AllK));                                       % Construct Space

% Risk-Neutral Density (Empirical)
Index_EMP = find((Smooth_AllK >= K_L0) & (Smooth_AllK <= K_R0));
Smooth_RND(Index_EMP) = Smooth_EMP_PDF((Smooth_K >= K_L0) & (Smooth_K <= K_R0));

% Generalized Extreme Value Function (GEV, Right Tail)
Index_GEV_R = find(Smooth_AllK >= K_R0);
try
    Smooth_RND(Index_GEV_R) = Smooth_GEV_R_PDF(Index_GEV_R);
catch
end

% Generalized Extreme Value Function (GEV, Left Tail)
Index_GEV_L = find(Smooth_AllK <= K_L0);
try
    Smooth_RND(Index_GEV_L) = Smooth_GEV_L_PDF(Index_GEV_L);
catch
end


%% Plot Figure: Risk-Neutral Density (Partial Density)

figure

plot(Smooth_AllK(Index_EMP), Smooth_RND(Index_EMP), ...
     'LineWidth', 2, ...
     'Color', 'r');     
hold on

set(gca,'XTick', min(Smooth_AllK):10:max(Smooth_AllK));
xlim([550 700]);
grid on

h1 = title(['Partial Estimated Risk-Neutral Density Function (Date = ' datestr(datenum(num2str(Target_Date), 'yyyymmdd'), 'mmm dd, yyyy') ', Time to Maturity = ' num2str(Target_TTM) ' Days)']);
h2 = legend('Empirical RND');
h3 = xlabel('Level of S&P 500 Index');
h4 = ylabel('Density');

set([h1 h2 h3 h4], ...
    'FontSize', 14, ...
    'FontWeight', 'bold', ...
    'LineWidth', 2)
clear h1 h2 h3 h4


%% Plot Figure: Risk-Neutral Density (Full Density)

figure

plot(Smooth_AllK(Index_EMP), Smooth_RND(Index_EMP), ...
     'LineWidth', 2, ...
     'Color', 'r');     
hold on

try
    plot(Smooth_AllK(Index_GEV_R), Smooth_RND(Index_GEV_R), ...
         'LineWidth', 2, ...
         'LineStyle', '--', ...
         'Color', 'b');     
    hold on
catch
end

try
    plot(Smooth_AllK(Index_GEV_L), Smooth_RND(Index_GEV_L), ...
         'LineWidth', 2, ...
         'LineStyle', ':', ...
         'Color', 'g');     
catch
end

set(gca,'XTick', min(Smooth_AllK):100:max(Smooth_AllK));
xlim([min(Smooth_AllK) max(Smooth_AllK)]);
grid on

h1 = title(['Full Estimated Risk-Neutral Density Function (Date = ' datestr(datenum(num2str(Target_Date), 'yyyymmdd'), 'mmm dd, yyyy') ', Time to Maturity = ' num2str(Target_TTM) ' Days)']);
h2 = legend('Empirical RND', 'Right GEV Tail', 'Left GEV Tail');
h3 = xlabel('Level of S&P 500 Index');
h4 = ylabel('Density');

set([h1 h2 h3 h4], ...
    'FontSize', 14, ...
    'FontWeight', 'bold', ...
    'LineWidth', 2)
clear h1 h2 h3 h4