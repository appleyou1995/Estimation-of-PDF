clear; clc


%% Setting

currentFolder = pwd;

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
Index = (Data(:, Index_Date)==Target_Date) & ...
        (Data(:, Index_TTM) >= Target_TTM);
Data = Data(Index, :);                                                     % Update: Data
clear Index

% Specific Time-to-Maturity 
Target_TTM = min(Data(:, Index_TTM));                                      % Update: Target_TTM

Index = Data(:, Index_TTM)==Target_TTM;
Data = Data(Index, :);                                                     % Update: Data
clear Index


%% Correction of Expiration Date (Time-to-Maturity)

% Saturday to Friday
dateTime = datetime(num2str(Data(:, Index_Date)), 'InputFormat', 'yyyyMMdd');
Date_EXP_WeekDay = weekday(dateTime + days(Data(:, Index_TTM)));
Data(Date_EXP_WeekDay==7, Index_TTM) = Data(Date_EXP_WeekDay==7, Index_TTM) - 1;   % Update: Data
clear Date_EXP_WeekDay       
    
% AM Settlement
% Calculation of Risk-Free Rate and Implied Volatility (IvyDB Reference Manual) 
Data(:, Index_TTM) = Data(:, Index_TTM) - 1;                               % Update: Data   

Target_TTM = Data(1, Index_TTM);


%% Data Filtering

% Bid >= 0.375
Index = Data(:, Index_OP_Bid) >= 0.375;
Data = Data(Index, :);

clear Index


% Out-of-the-Money (OTM) and Around At-the-Money (ATM) Options
Index = (Data(:, Index_CPFlag) == 1) & ...
        (Data(:, Index_K) >= Data(:, Index_S)) & ...
        (Data(:, Index_K) <= Data(:, Index_S));
Data(Index, Index_CPFlag) = 31;                                            % Update: Data (CP Flag)           
clear Index

Index = (Data(:, Index_CPFlag)==2) & ...
        (Data(:, Index_K) >= (Data(:, Index_S))) & ...
        (Data(:, Index_K) <= (Data(:, Index_S)));
Data(Index, Index_CPFlag) = 32;                                            % Update: Data (CP Flag)           
clear Index
          
Index = ((Data(:, Index_CPFlag)==1) & (Data(:, Index_K) >= Data(:, Index_S))) | ...
        ((Data(:, Index_CPFlag)==2) & (Data(:, Index_K) <= Data(:, Index_S))) | ...
        (fix(Data(:, Index_CPFlag) / 10)==3);      
Data = Data(Index, :);                                                     % Update: Data 
clear Index 


%% Calculation of Dividend Yield

% [1. SecID | 2. Date (YYYYMMDD) | 3. Dividend Yield (Annualized)]

FileName = 'IndexDivYield19962019.txt';
Data_DY = load(['D:\Google\我的雲端硬碟\學術｜研究與論文\論文著作\CDI Method\Data\99 姿穎學姊提供\20240417\' FileName]);
clear FileName

Data(:, Index_DY) = Data_DY(Data_DY(:, Index_Date)==Target_Date, end);     % Annualized
% clear Data_DY


%% Calculation of Risk-Free Rate

% [1. Date (YYYYMMDD) | 2. TTM (Days) | 3. Risk-Free Rate (Annualized)]

FileName = 'RiskFreeRate19962019.txt';
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
[Call, Put] = blsprice(S0, Smooth_K, RF, TTM, Smooth_IV, DY); 

% clear Smooth_IV
% clear S0 DY


%% Define Color (LaTeX Beamer Theme - Metropolis)

mRed        = '#e74c3c';
mDarkRed    = '#b22222';
mLightBlue  = '#3279a8';
mDarkBlue   = '#2c3e50';
mDarkGreen  = '#4b8b3b';
mOrange     = '#f39c12';
mBackground = '#FAFAFA';


%% Plot Figure: Smoothed Implied Volatility (Fourth Order Spline Approximation)

ImpliedVolatilities = figure;

set(gcf, 'Color', mBackground);

plot(Smooth_K, Smooth_IV, ...
     'LineWidth', 2, ...
     'Color', mOrange);   
hold on

plot(Data(:, Index_K), Data(:, Index_IV), ...
     'LineStyle', 'none', ...     
     'LineWidth', 2, ...
     'Marker', 'o', ...     
     'MarkerSize', 6, ...
     'Color', mDarkBlue); 
grid on

h1 = title('Implied Volatilities with Interpolation and Smoothing');
h2 = legend('4th Degree Spline Interpolation (Without Knot)', 'Raw Option Data');
h3 = xlabel('Strike Price');
h4 = ylabel('Implied Volatility');

% Adjust the title position
titlePosition = get(h1, 'Position');
titlePosition(2) = titlePosition(2) + 0.005; % Adjust this value as needed
set(h1, 'Position', titlePosition);

set([h1 h2 h3 h4], ...
    'FontName', 'Fira Sans', ...
    'FontSize', 10, ...
    'LineWidth', 1)

% Customize axes
ax = gca;
set(ax, 'FontName', 'Fira Sans Light', 'FontSize', 12);

set(gcf,'InvertHardcopy','off')
clear h1 h2 h3 h4

% Save the figure
saveas(ImpliedVolatilities, fullfile(currentFolder, 'Figure', 'ImpliedVolatilities.png'));
saveas(ImpliedVolatilities, fullfile(currentFolder, 'Figure', 'ImpliedVolatilities.fig'));


%% Plot Figure: Smoothed Option Price

Option_Price = figure;

set(gcf, 'Color', mBackground);

% Customize axes
ax = gca;
set(ax, 'FontName', 'Fira Sans Light', 'FontSize', 12);

S0_line = Data(1, Index_S);
line([S0_line S0_line], [0.2 41], 'Color', mDarkRed, 'LineWidth', 1);
hold on

text(S0_line, 42, sprintf('%.2f', S0_line), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 11, 'FontName', 'Fira Sans Medium', 'Color', mDarkRed);

plot(Smooth_K, Call, ...
     'LineWidth', 2, ...
     'LineStyle', '--', ...
     'Color', mLightBlue);   

plot(Smooth_K, Put, ...
     'LineWidth', 2, ...
     'LineStyle', ':', ...
     'Color', mDarkGreen);  

plot(Data(:, Index_K), Data(:, Index_OP_Bid), ...
     'LineStyle', 'none', ...     
     'LineWidth', 2, ...
     'Marker', 'o', ...     
     'MarkerSize', 6, ...
     'Color', mDarkBlue); 

grid on

h1 = title('Option Prices via Implied Volatility Smoothing');
h2 = legend('Spot Price', 'Call Option Price', 'Put Option Price', 'Raw Option Price Data');
h3 = xlabel('Strike Price');
h4 = ylabel('Option Price');

% Adjust the title position
titlePosition = get(h1, 'Position');
titlePosition(2) = titlePosition(2) + 0.05; % Adjust this value as needed
set(h1, 'Position', titlePosition);

set([h1 h2 h3 h4], ...
    'FontName', 'Fira Sans', ...
    'FontSize', 10, ...
    'LineWidth', 1)

ax = gca;
set(ax, 'FontName', 'Fira Sans Light', 'FontSize', 12, 'box', 'on');

set(gcf,'InvertHardcopy','off')
clear h1 h2 h3 h4

% Save the figure
saveas(Option_Price, fullfile(currentFolder, 'Figure', 'Option_Price.png'));
saveas(Option_Price, fullfile(currentFolder, 'Figure', 'Option_Price.fig'));


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


%% Plot Figure: Risk-Neutral Density (Partial Density)

Q_measure_PDF_Partial = figure;

set(gcf, 'Color', mBackground);

plot(Smooth_K, Smooth_EMP_PDF, ...
     'LineWidth', 2, ...
     'Color', mDarkRed);     
hold on

set(gca,'XTick', min(Smooth_AllK):20:max(Smooth_AllK));
xlim([540 690]);
xtickformat('%.0f');
grid on

h1 = title('Partial Estimated Risk-Neutral Density Function');
h2 = legend('Empirical RND', 'Location', 'northwest');
h3 = xlabel('Level of S&P 500 Index');
h4 = ylabel('Probability Density');

% Adjust the title position
titlePosition = get(h1, 'Position');
titlePosition(2) = titlePosition(2) + 0.0005; % Adjust this value as needed
set(h1, 'Position', titlePosition);

set([h1 h2 h3 h4], ...
    'FontName', 'Fira Sans', ...
    'FontSize', 10, ...
    'LineWidth', 2)

% Customize axes
ax = gca;
set(ax, 'FontName', 'Fira Sans Light', 'FontSize', 12);

set(gcf,'InvertHardcopy','off')

clear h1 h2 h3 h4

% Save the figure
saveas(Q_measure_PDF_Partial, fullfile(currentFolder, 'Figure', 'Q_measure_PDF_Partial.png'));
saveas(Q_measure_PDF_Partial, fullfile(currentFolder, 'Figure', 'Q_measure_PDF_Partial.fig')); 


%% Plot Figure: CDF

Q_measure_CDF = figure;

set(gcf, 'Color', mBackground);

plot(Smooth_K, Smooth_EMP_CDF, ...
     'LineWidth', 2, ...
     'Color', mDarkRed);  
hold on

set(gca,'XTick', min(Smooth_AllK):20:max(Smooth_AllK));
xlim([540 690]);
xtickformat('%.0f');
grid on

h1 = title('Empirical Cumulative Distribution Function');
h2 = legend('Empirical CDF', 'Location', 'northwest');
h3 = xlabel('Level of S&P 500 Index');
h4 = ylabel('Cummulated Density');

% Adjust the title position
titlePosition = get(h1, 'Position');
titlePosition(2) = titlePosition(2) + 0.015; % Adjust this value as needed
set(h1, 'Position', titlePosition);

set([h1 h2 h3 h4], ...
    'FontName', 'Fira Sans', ...
    'FontSize', 10, ...
    'LineWidth', 2)

% Customize axes
ax = gca;
set(ax, 'FontName', 'Fira Sans Light', 'FontSize', 12);

set(gcf,'InvertHardcopy','off')

clear h1 h2 h3 h4

% Save the figure
saveas(Q_measure_CDF, fullfile(currentFolder, 'Figure', 'Q_measure_CDF.png'));
saveas(Q_measure_CDF, fullfile(currentFolder, 'Figure', 'Q_measure_CDF.fig'));


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
BP_R0 = BP_R1 - 0.03;                                                      % Update: BP_R0 
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

% clear EMP_CDF_R0 EMP_PDF_R0 EMP_PDF_R1
% clear mu sigma k sol


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

% clear EMP_CDF_L0 EMP_PDF_L0 EMP_PDF_L1  
% clear mu sigma k sol


%% Plot Figure: Right Tail

GEV_Right_Tail = figure;

set(gcf, 'Color', mBackground);

% Risk-Neutral Density (Empirical)
plot(Smooth_K, Smooth_EMP_PDF, ...
     'LineWidth', 2, ...
     'Color', mDarkRed);     
hold on
 
% Generalized Extreme Value Function (GEV, Right Tail)
try
    plot(Smooth_AllK, Smooth_GEV_R_PDF, ...
         'LineWidth', 2, ...
         'LineStyle', '--', ...
         'Color', mLightBlue);    
    hold on
catch
end

set(gca,'XTick', min(Smooth_AllK):50:max(Smooth_AllK));
xlim([450 900]);
ylim([0 0.025]);
xtickformat('%.0f');
grid on

h1 = title('Risk-Neutral Density and Fitted GEV Right Tail Functions');
h2 = legend('Empirical RND', 'Right Tail GEV Function');
h3 = xlabel('Level of S&P 500 Index');
h4 = ylabel('Density');

% Adjust the title position
titlePosition = get(h1, 'Position');
titlePosition(2) = titlePosition(2) + 0.0005;
set(h1, 'Position', titlePosition);

set([h1 h2 h3 h4], ...
    'FontName', 'Fira Sans', ...
    'FontSize', 10, ...
    'LineWidth', 1)

% Customize axes
ax = gca;
set(ax, 'FontName', 'Fira Sans Light', 'FontSize', 12);

set(gcf,'InvertHardcopy','off')

clear h1 h2 h3 h4

% Save the figure
saveas(GEV_Right_Tail, fullfile(currentFolder, 'Figure', 'GEV_Right_Tail.png'));
saveas(GEV_Right_Tail, fullfile(currentFolder, 'Figure', 'GEV_Right_Tail.fig'));


%% Plot Figure: Left Tail

GEV_Left_Tail = figure;

set(gcf, 'Color', mBackground);

% Risk-Neutral Density (Empirical)
plot(Smooth_K, Smooth_EMP_PDF, ...
     'LineWidth', 2, ...
     'Color', mDarkRed);     
hold on 

% Generalized Extreme Value Function (GEV, Left Tail)
try
    plot(Smooth_AllK, Smooth_GEV_L_PDF, ...
         'LineWidth', 2, ...
         'LineStyle', ':', ...
         'Color', mDarkGreen);    
    hold on
catch
end 

set(gca,'XTick', min(Smooth_AllK):50:max(Smooth_AllK));
xlim([450 900]);
ylim([0 0.025]);
xtickformat('%.0f');
grid on

h1 = title('Risk-Neutral Density and Fitted GEV Left Tail Functions');
h2 = legend('Empirical RND', 'Left Tail GEV Function');
h3 = xlabel('Level of S&P 500 Index');
h4 = ylabel('Density');

% Adjust the title position
titlePosition = get(h1, 'Position');
titlePosition(2) = titlePosition(2) + 0.0005;
set(h1, 'Position', titlePosition);

set([h1 h2 h3 h4], ...
    'FontName', 'Fira Sans', ...
    'FontSize', 10, ...
    'LineWidth', 1)

% Customize axes
ax = gca;
set(ax, 'FontName', 'Fira Sans Light', 'FontSize', 12);

set(gcf,'InvertHardcopy','off')

clear h1 h2 h3 h4

% Save the figure
saveas(GEV_Left_Tail, fullfile(currentFolder, 'Figure', 'GEV_Left_Tail.png'));
saveas(GEV_Left_Tail, fullfile(currentFolder, 'Figure', 'GEV_Left_Tail.fig'));


%% Plot Figure

figure

% Risk-Neutral Density (Empirical)
plot(Smooth_K, Smooth_EMP_PDF, ...
     'LineWidth', 2, ...
     'Color', mDarkRed);     
hold on
 
% Generalized Extreme Value Function (GEV, Right Tail)
try
    plot(Smooth_AllK, Smooth_GEV_R_PDF, ...
         'LineWidth', 2, ...
         'LineStyle', '--', ...
         'Color', mLightBlue);    
    hold on
catch
end  

% Generalized Extreme Value Function (GEV, Left Tail)
try
    plot(Smooth_AllK, Smooth_GEV_L_PDF, ...
         'LineWidth', 2, ...
         'LineStyle', ':', ...
         'Color', mDarkGreen);    
    hold on
catch
end

% Connection Points
plot(Smooth_K(Smooth_K==K_R0), Smooth_EMP_PDF(Smooth_K==K_R0), ...
     'LineWidth', 2, ...
     'Marker', 'o', ...
     'Color', mOrange);    
text(Smooth_K(Smooth_K==K_R0), Smooth_EMP_PDF(Smooth_K==K_R0), ...
     ['  \leftarrow ' num2str(100 * roundn(BP_R0, - 4)) '%'], ...
     'HorizontalAlignment', 'left', ...
     'VerticalAlignment', 'baseline', ...
     'FontSize', 11, ...
     'FontName', 'Fira Sans Medium', ...
     'LineWidth', 2, ...
     'Color', mOrange);        
hold on

plot(Smooth_K(Smooth_K==K_R1), Smooth_EMP_PDF(Smooth_K==K_R1), ...
     'LineWidth', 2, ...
     'Marker', 'o', ...
     'Color', mOrange);    
text(Smooth_K(Smooth_K==K_R1), Smooth_EMP_PDF(Smooth_K==K_R1), ...
     ['  \leftarrow ' num2str(100 * roundn(BP_R1, - 4)) '%'], ...
     'HorizontalAlignment', 'left', ...
     'VerticalAlignment', 'baseline', ...
     'FontSize', 11, ...
     'FontName', 'Fira Sans Medium', ...
     'LineWidth', 2, ...
     'Color', mOrange);       
hold on

plot(Smooth_K(Smooth_K==K_L0), Smooth_EMP_PDF(Smooth_K==K_L0), ...
     'LineWidth', 2, ...
     'Marker', 'o', ...
     'Color', mOrange);    
text(Smooth_K(Smooth_K==K_L0), Smooth_EMP_PDF(Smooth_K==K_L0), ...
     [num2str(100 * roundn(BP_L0, - 4)) '% \rightarrow  '], ...
     'HorizontalAlignment', 'right', ...
     'VerticalAlignment', 'baseline', ...
     'FontSize', 11, ...
     'FontName', 'Fira Sans Medium', ...
     'LineWidth', 2, ...
     'Color', mOrange);     
hold on

plot(Smooth_K(Smooth_K==K_L1), Smooth_EMP_PDF(Smooth_K==K_L1), ...
     'LineWidth', 2, ...
     'Marker', 'o', ...
     'Color', mOrange);  
text(Smooth_K(Smooth_K==K_L1), Smooth_EMP_PDF(Smooth_K==K_L1), ...
     [num2str(100 * roundn(BP_L1, - 4)) '% \rightarrow  '], ...
     'HorizontalAlignment', 'right', ...
     'VerticalAlignment', 'baseline', ...
     'FontSize', 11, ...
     'FontName', 'Fira Sans Medium', ...
     'LineWidth', 2, ...
     'Color', mOrange);   

set(gca,'XTick', min(Smooth_AllK):100:max(Smooth_AllK));
xlim([min(Smooth_AllK) max(Smooth_AllK)]);
xtickformat('%.0f');
grid on

h1 = title('Risk-Neutral Density and Fitted GEV Tail Functions');
h2 = legend('Empirical RND', 'Right Tail GEV Function', 'Left Tail GEV Function', 'Connection Points');
h3 = xlabel('Level of S&P 500 Index');
h4 = ylabel('Density');

set([h1 h2 h3 h4], ...
    'FontName', 'Fira Sans', ...
    'FontSize', 10, ...
    'LineWidth', 1)

% Customize axes
ax = gca;
set(ax, 'FontName', 'Fira Sans', 'FontSize', 12);

set(gcf,'InvertHardcopy','off')

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


%% Plot Figure: Risk-Neutral Density (Full Density)

Q_measure_PDF_Full = figure;

set(gcf, 'Color', mBackground);

plot(Smooth_AllK(Index_EMP), Smooth_RND(Index_EMP), ...
     'LineWidth', 2, ...
     'Color', mDarkRed);     
hold on

try
    plot(Smooth_AllK(Index_GEV_R), Smooth_RND(Index_GEV_R), ...
         'LineWidth', 2, ...
         'LineStyle', '--', ...
         'Color', mLightBlue);     
    hold on
catch
end

try
    plot(Smooth_AllK(Index_GEV_L), Smooth_RND(Index_GEV_L), ...
         'LineWidth', 2, ...
         'LineStyle', ':', ...
         'Color', mDarkGreen);     
catch
end

set(gca,'XTick', min(Smooth_AllK):100:max(Smooth_AllK));
xlim([min(Smooth_AllK) max(Smooth_AllK)]);
xtickformat('%.0f');
grid on

h1 = title('Full Estimated Risk-Neutral Density Function');
h2 = legend('Empirical RND', 'Right GEV Tail', 'Left GEV Tail');
h3 = xlabel('Level of S&P 500 Index');
h4 = ylabel('Probability Density');

% Adjust the title position
titlePosition = get(h1, 'Position');
titlePosition(2) = titlePosition(2) + 0.0005;
set(h1, 'Position', titlePosition);

set([h1 h2 h3 h4], ...
    'FontName', 'Fira Sans', ...
    'FontSize', 10, ...
    'LineWidth', 2)

% Customize axes
ax = gca;
set(ax, 'FontName', 'Fira Sans Light', 'FontSize', 12);

set(gcf,'InvertHardcopy','off')

clear h1 h2 h3 h4

% Save the figure
saveas(Q_measure_PDF_Full, fullfile(currentFolder, 'Figure', 'Q_measure_PDF_Full.png'));
saveas(Q_measure_PDF_Full, fullfile(currentFolder, 'Figure', 'Q_measure_PDF_Full.fig'));


%% [CARA] Constant Absolute Risk Aversion Utility Function

alpha_1 = 0.01;
alpha_2 = 0.02;
alpha_3 = 0.03;
alpha_4 = 0.04;

% [Negative Exponential Utility] Setting Numerator with Different Alpha

utility_alpha_1 = (1 / alpha_1) .* exp(alpha_1 .* Smooth_AllK) .* Smooth_RND;
utility_alpha_2 = (1 / alpha_2) .* exp(alpha_2 .* Smooth_AllK) .* Smooth_RND;
utility_alpha_3 = (1 / alpha_3) .* exp(alpha_3 .* Smooth_AllK) .* Smooth_RND;
utility_alpha_4 = (1 / alpha_4) .* exp(alpha_4 .* Smooth_AllK) .* Smooth_RND;

% [Negative Exponential Utility] Real-World Densities
P_Density_alpha_1 = utility_alpha_1 ./ trapz(Smooth_AllK, utility_alpha_1);
P_Density_alpha_2 = utility_alpha_2 ./ trapz(Smooth_AllK, utility_alpha_2);
P_Density_alpha_3 = utility_alpha_3 ./ trapz(Smooth_AllK, utility_alpha_3);
P_Density_alpha_4 = utility_alpha_4 ./ trapz(Smooth_AllK, utility_alpha_4);


%% [CARA] Plot Figure: Physical Density (Full Density)

P_measure_PDF_CARA_Full = figure;

set(gcf, 'Color', mBackground);

plot(Smooth_AllK, P_Density_alpha_1, '-', 'LineWidth', 2);
hold on;
plot(Smooth_AllK, P_Density_alpha_2, '-', 'LineWidth', 2);
plot(Smooth_AllK, P_Density_alpha_3, '-', 'LineWidth', 2);
plot(Smooth_AllK, P_Density_alpha_4, '-', 'LineWidth', 2);
hold off;

set(gca,'XTick', min(Smooth_AllK):100:max(Smooth_AllK));
set(gca,'YTick', 0:0.002:0.02);
xlim([min(Smooth_AllK) max(Smooth_AllK)]);
xtickformat('%.0f');
grid on

legendText1 = sprintf('alpha=%.2f', alpha_1);
legendText2 = sprintf('alpha=%.2f', alpha_2);
legendText3 = sprintf('alpha=%.2f', alpha_3);
legendText4 = sprintf('alpha=%.2f', alpha_4);

h1 = title('Full Physical Density Function (Negative Exponential Utility)');
h2 = legend(legendText1, legendText2, legendText3, legendText4);
h3 = xlabel('Level of S&P 500 Index');
h4 = ylabel('Probability Density');

% Adjust the title position
titlePosition = get(h1, 'Position');
titlePosition(2) = titlePosition(2) + 0.0005;
set(h1, 'Position', titlePosition);

set([h1 h2 h3 h4], ...
    'FontName', 'Fira Sans', ...
    'FontSize', 10, ...
    'LineWidth', 2)

% Customize axes
ax = gca;
set(ax, 'FontName', 'Fira Sans Light', 'FontSize', 12);

set(gcf,'InvertHardcopy','off')

clear h1 h2 h3 h4 legendText1 legendText2 legendText3 legendText4

% Save the figure
saveas(P_measure_PDF_CARA_Full, fullfile(currentFolder, 'Figure', 'P_measure_PDF_CARA_Full.png'));
saveas(P_measure_PDF_CARA_Full, fullfile(currentFolder, 'Figure', 'P_measure_PDF_CARA_Full.fig'));


%% [CARA] Plot Figure: Physical Density (Partial Density)

P_measure_PDF_CARA_Partial = figure;

set(gcf, 'Color', mBackground);

plot(Smooth_AllK, P_Density_alpha_1, '-', 'LineWidth', 2);
hold on;
plot(Smooth_AllK, P_Density_alpha_2, '-', 'LineWidth', 2);
plot(Smooth_AllK, P_Density_alpha_3, '-', 'LineWidth', 2);
plot(Smooth_AllK, P_Density_alpha_4, '-', 'LineWidth', 2);
hold off;

set(gca,'XTick', min(Smooth_AllK):40:max(Smooth_AllK));
set(gca,'YTick', 0:0.002:0.02);
xlim([600, 740]);
ylim([0.01, 0.02]);
xtickformat('%.0f');
grid on

legendText1 = sprintf('alpha=%.2f', alpha_1);
legendText2 = sprintf('alpha=%.2f', alpha_2);
legendText3 = sprintf('alpha=%.2f', alpha_3);
legendText4 = sprintf('alpha=%.2f', alpha_4);

h1 = title('Partial Physical Density Function (Negative Exponential Utility)');
h2 = legend(legendText1, legendText2, legendText3, legendText4);
h3 = xlabel('Level of S&P 500 Index');
h4 = ylabel('Probability Density');

% Adjust the title position
titlePosition = get(h1, 'Position');
titlePosition(2) = titlePosition(2) + 0.0002;
set(h1, 'Position', titlePosition);

set([h1 h2 h3 h4], ...
    'FontName', 'Fira Sans', ...
    'FontSize', 10, ...
    'LineWidth', 2)

% Customize axes
ax = gca;
set(ax, 'FontName', 'Fira Sans Light', 'FontSize', 12);

set(gcf,'InvertHardcopy','off')

clear h1 h2 h3 h4 legendText1 legendText2 legendText3 legendText4

% Save the figure
saveas(P_measure_PDF_CARA_Partial, fullfile(currentFolder, 'Figure', 'P_measure_PDF_CARA_Partial.png'));
saveas(P_measure_PDF_CARA_Partial, fullfile(currentFolder, 'Figure', 'P_measure_PDF_CARA_Partial.fig'));


%% [CARA] Plot Figure: Physical Density (Zoom-in Density)

P_measure_PDF_CARA_Zoom = figure;

set(gcf, 'Color', mBackground);

plot(Smooth_AllK, P_Density_alpha_1, '-', 'LineWidth', 2, 'Color', mDarkGreen);
hold on;
plot(Smooth_AllK, P_Density_alpha_2, '-', 'LineWidth', 2, 'Color', mOrange);
plot(Smooth_AllK, P_Density_alpha_3, '-', 'LineWidth', 2, 'Color', mLightBlue);
plot(Smooth_AllK, P_Density_alpha_4, '-', 'LineWidth', 2, 'Color', mDarkRed);
hold off;

set(gca,'XTick', min(Smooth_AllK):100:max(Smooth_AllK));
set(gca,'YTick', 0:0.002:0.02);
xlim([min(Smooth_AllK) max(Smooth_AllK)]);
xtickformat('%.0f');
grid on

legendText1 = sprintf('alpha=%.2f', alpha_1);
legendText2 = sprintf('alpha=%.2f', alpha_2);
legendText3 = sprintf('alpha=%.2f', alpha_3);
legendText4 = sprintf('alpha=%.2f', alpha_4);

h1 = title('Full Physical Density Function (Negative Exponential Utility)');
h2 = legend(legendText1, legendText2, legendText3, legendText4);
h3 = xlabel('Level of S&P 500 Index');
h4 = ylabel('Probability Density');

% Adjust the title position
titlePosition = get(h1, 'Position');
titlePosition(2) = titlePosition(2) + 0.0005;
set(h1, 'Position', titlePosition);

set([h1 h2 h3 h4], ...
    'FontName', 'Fira Sans', ...
    'FontSize', 10, ...
    'LineWidth', 2)

% Customize axes
ax = gca;
set(ax, 'FontName', 'Fira Sans Light', 'FontSize', 12);

% Add inset axes for zoomed-in plot
inset_axes = axes('Position',[0.55 0.3 0.3 0.35]);
box on
plot(Smooth_AllK, P_Density_alpha_1, '-', 'LineWidth', 2, 'Color', mDarkGreen);
hold on;
plot(Smooth_AllK, P_Density_alpha_2, '-', 'LineWidth', 2, 'Color', mOrange);
plot(Smooth_AllK, P_Density_alpha_3, '-', 'LineWidth', 2, 'Color', mLightBlue);
plot(Smooth_AllK, P_Density_alpha_4, '-', 'LineWidth', 2, 'Color', mDarkRed);
hold off;
set(gca, 'XTick', 620:20:700, 'YTick', 0.01:0.002:0.02);
xlim([620, 700]);
ylim([0.01, 0.02]);
grid on

set(inset_axes, 'FontName', 'Fira Sans Light', 'FontSize', 8);

set(gcf,'InvertHardcopy','off')

clear h1 h2 h3 h4 legendText1 legendText2 legendText3 legendText4

% Save the figure
saveas(P_measure_PDF_CARA_Zoom, fullfile(currentFolder, 'Figure', 'P_measure_PDF_CARA_Zoom.png'));
saveas(P_measure_PDF_CARA_Zoom, fullfile(currentFolder, 'Figure', 'P_measure_PDF_CARA_Zoom.fig'));


%% [CRRA] Constant Relative Risk Aversion Utility Function

gamma_1 = 1;
gamma_2 = 2;
gamma_3 = 3;
gamma_4 = 4;

% [Power Utility] Setting Numerator with Different Gamma
utility_gamma_1 = Smooth_AllK .^ gamma_1 .* Smooth_RND;
utility_gamma_2 = Smooth_AllK .^ gamma_2 .* Smooth_RND;
utility_gamma_3 = Smooth_AllK .^ gamma_3 .* Smooth_RND;
utility_gamma_4 = Smooth_AllK .^ gamma_4 .* Smooth_RND;

% [Power Utility] Real-World Densities
P_Density_gamma_1 = utility_gamma_1 ./ trapz(Smooth_AllK, utility_gamma_1);
P_Density_gamma_2 = utility_gamma_2 ./ trapz(Smooth_AllK, utility_gamma_2);
P_Density_gamma_3 = utility_gamma_3 ./ trapz(Smooth_AllK, utility_gamma_3);
P_Density_gamma_4 = utility_gamma_4 ./ trapz(Smooth_AllK, utility_gamma_4);


%% [CRRA] Plot Figure: Physical Density (Full Density)

P_measure_PDF_CRRA_Full = figure;

set(gcf, 'Color', mBackground);

plot(Smooth_AllK, P_Density_gamma_1, '-', 'LineWidth', 2);
hold on;
plot(Smooth_AllK, P_Density_gamma_2, '-', 'LineWidth', 2);
plot(Smooth_AllK, P_Density_gamma_3, '-', 'LineWidth', 2);
plot(Smooth_AllK, P_Density_gamma_4, '-', 'LineWidth', 2);
hold off;

set(gca,'XTick', min(Smooth_AllK):100:max(Smooth_AllK));
set(gca,'YTick', 0:0.002:0.02);
xlim([min(Smooth_AllK) max(Smooth_AllK)]);
xtickformat('%.0f');
grid on

legendText1 = sprintf('gamma=%.0f', gamma_1);
legendText2 = sprintf('gamma=%.0f', gamma_2);
legendText3 = sprintf('gamma=%.0f', gamma_3);
legendText4 = sprintf('gamma=%.0f', gamma_4);

h1 = title('Full Physical Density Function (Power Utility)');
h2 = legend(legendText1, legendText2, legendText3, legendText4);
h3 = xlabel('Level of S&P 500 Index');
h4 = ylabel('Probability Density');

% Adjust the title position
titlePosition = get(h1, 'Position');
titlePosition(2) = titlePosition(2) + 0.0005;
set(h1, 'Position', titlePosition);

set([h1 h2 h3 h4], ...
    'FontName', 'Fira Sans', ...
    'FontSize', 10, ...
    'LineWidth', 2)

% Customize axes
ax = gca;
set(ax, 'FontName', 'Fira Sans Light', 'FontSize', 12);

set(gcf,'InvertHardcopy','off')

clear h1 h2 h3 h4 legendText1 legendText2 legendText3 legendText4

% Save the figure
saveas(P_measure_PDF_CRRA_Full, fullfile(currentFolder, 'Figure', 'P_measure_PDF_CRRA_Full.png'));
saveas(P_measure_PDF_CRRA_Full, fullfile(currentFolder, 'Figure', 'P_measure_PDF_CRRA_Full.fig'));


%% [CRRA] Plot Figure: Physical Density (Partial Density)

P_measure_PDF_CRRA_Partial = figure;

set(gcf, 'Color', mBackground);

plot(Smooth_AllK, P_Density_gamma_1, '-', 'LineWidth', 2);
hold on;
plot(Smooth_AllK, P_Density_gamma_2, '-', 'LineWidth', 2);
plot(Smooth_AllK, P_Density_gamma_3, '-', 'LineWidth', 2);
plot(Smooth_AllK, P_Density_gamma_4, '-', 'LineWidth', 2);
hold off;

set(gca,'XTick', min(Smooth_AllK):10:max(Smooth_AllK));
set(gca,'YTick', 0:0.001:0.02);
xlim([600, 700]);
ylim([0.014, 0.017]);
xtickformat('%.0f');
grid on

legendText1 = sprintf('gamma=%.0f', gamma_1);
legendText2 = sprintf('gamma=%.0f', gamma_2);
legendText3 = sprintf('gamma=%.0f', gamma_3);
legendText4 = sprintf('gamma=%.0f', gamma_4);

h1 = title('Partial Physical Density Function (Power Utility)');
h2 = legend(legendText1, legendText2, legendText3, legendText4);
h3 = xlabel('Level of S&P 500 Index');
h4 = ylabel('Probability Density');

% Adjust the title position
titlePosition = get(h1, 'Position');
titlePosition(2) = titlePosition(2) + 0.0001;
set(h1, 'Position', titlePosition);

set([h1 h2 h3 h4], ...
    'FontName', 'Fira Sans', ...
    'FontSize', 10, ...
    'LineWidth', 2)

% Customize axes
ax = gca;
set(ax, 'FontName', 'Fira Sans Light', 'FontSize', 12);

set(gcf,'InvertHardcopy','off')

clear h1 h2 h3 h4 legendText1 legendText2 legendText3 legendText4

% Save the figure
saveas(P_measure_PDF_CRRA_Partial, fullfile(currentFolder, 'Figure', 'P_measure_PDF_CRRA_Partial.png'));
saveas(P_measure_PDF_CRRA_Partial, fullfile(currentFolder, 'Figure', 'P_measure_PDF_CRRA_Partial.fig'));


%% [CRRA] Plot Figure: Physical Density (Zoom-in Density)

P_measure_PDF_CRRA_Zoom = figure;

set(gcf, 'Color', mBackground);

plot(Smooth_AllK, P_Density_gamma_1, '-', 'LineWidth', 2, 'Color', mDarkGreen);
hold on;
plot(Smooth_AllK, P_Density_gamma_2, '-', 'LineWidth', 2, 'Color', mOrange);
plot(Smooth_AllK, P_Density_gamma_3, '-', 'LineWidth', 2, 'Color', mLightBlue);
plot(Smooth_AllK, P_Density_gamma_4, '-', 'LineWidth', 2, 'Color', mDarkRed);
hold off;

set(gca,'XTick', min(Smooth_AllK):100:max(Smooth_AllK));
set(gca,'YTick', 0:0.002:0.02);
xlim([min(Smooth_AllK) max(Smooth_AllK)]);
xtickformat('%.0f');
grid on

legendText1 = sprintf('gamma=%.0f', gamma_1);
legendText2 = sprintf('gamma=%.0f', gamma_2);
legendText3 = sprintf('gamma=%.0f', gamma_3);
legendText4 = sprintf('gamma=%.0f', gamma_4);

h1 = title('Full Physical Density Function (Power Utility)');
h2 = legend(legendText1, legendText2, legendText3, legendText4);
h3 = xlabel('Level of S&P 500 Index');
h4 = ylabel('Probability Density');

% Adjust the title position
titlePosition = get(h1, 'Position');
titlePosition(2) = titlePosition(2) + 0.0005;
set(h1, 'Position', titlePosition);

set([h1 h2 h3 h4], ...
    'FontName', 'Fira Sans', ...
    'FontSize', 10, ...
    'LineWidth', 2)

% Customize axes
ax = gca;
set(ax, 'FontName', 'Fira Sans Light', 'FontSize', 12);

% Add inset axes for zoomed-in plot
inset_axes = axes('Position',[0.55 0.3 0.3 0.35]);
box on
plot(Smooth_AllK, P_Density_gamma_1, '-', 'LineWidth', 2, 'Color', mDarkGreen);
hold on;
plot(Smooth_AllK, P_Density_gamma_2, '-', 'LineWidth', 2, 'Color', mOrange);
plot(Smooth_AllK, P_Density_gamma_3, '-', 'LineWidth', 2, 'Color', mLightBlue);
plot(Smooth_AllK, P_Density_gamma_4, '-', 'LineWidth', 2, 'Color', mDarkRed);
hold off;
set(gca, 'XTick', 630:10:670, 'YTick', 0.014:0.001:0.017);
xlim([630, 670]);
ylim([0.014, 0.017]);
grid on

set(inset_axes, 'FontName', 'Fira Sans Light', 'FontSize', 8);

set(gcf,'InvertHardcopy','off')

clear h1 h2 h3 h4 legendText1 legendText2 legendText3 legendText4

% Save the figure
saveas(P_measure_PDF_CRRA_Zoom, fullfile(currentFolder, 'Figure', 'P_measure_PDF_CRRA_Zoom.png'));
saveas(P_measure_PDF_CRRA_Zoom, fullfile(currentFolder, 'Figure', 'P_measure_PDF_CRRA_Zoom.fig'));


%% Calculate Statistical Moments of Return

S_t = Data(1, Index_S);
S_T = Smooth_AllK;
R_T = S_T / S_t;

RND_mean     = trapz(Smooth_AllK, R_T .* Smooth_RND);
RND_variance = trapz(Smooth_AllK, ((R_T - RND_mean) .^ 2) .* Smooth_RND);
RND_std_dev  = sqrt(RND_variance);
RND_skewness = trapz(Smooth_AllK, (((R_T - RND_mean) / RND_std_dev) .^ 3) .* Smooth_RND);
RND_kurtosis = trapz(Smooth_AllK, (((R_T - RND_mean) / RND_std_dev) .^ 4) .* Smooth_RND) - 3;


%%  [Negative Exponential Utility]

% alpha = 0.01
alpha_1_mean     = trapz(Smooth_AllK, R_T .* P_Density_alpha_1);
alpha_1_variance = trapz(Smooth_AllK, ((R_T - alpha_1_mean) .^ 2) .* P_Density_alpha_1);
alpha_1_std_dev  = sqrt(alpha_1_variance);
alpha_1_skewness = trapz(Smooth_AllK, (((R_T - alpha_1_mean) / alpha_1_std_dev) .^ 3) .* P_Density_alpha_1);
alpha_1_kurtosis = trapz(Smooth_AllK, (((R_T - alpha_1_mean) / alpha_1_std_dev) .^ 4) .* P_Density_alpha_1) - 3;


% alpha = 0.02
alpha_2_mean     = trapz(Smooth_AllK, R_T .* P_Density_alpha_2);
alpha_2_variance = trapz(Smooth_AllK, ((R_T - alpha_2_mean) .^ 2) .* P_Density_alpha_2);
alpha_2_std_dev  = sqrt(alpha_2_variance);
alpha_2_skewness = trapz(Smooth_AllK, (((R_T - alpha_2_mean) / alpha_2_std_dev) .^ 3) .* P_Density_alpha_2);
alpha_2_kurtosis = trapz(Smooth_AllK, (((R_T - alpha_2_mean) / alpha_2_std_dev) .^ 4) .* P_Density_alpha_2) - 3;


% alpha = 0.03
alpha_3_mean     = trapz(Smooth_AllK, R_T .* P_Density_alpha_3);
alpha_3_variance = trapz(Smooth_AllK, ((R_T - alpha_3_mean) .^ 2) .* P_Density_alpha_3);
alpha_3_std_dev  = sqrt(alpha_3_variance);
alpha_3_skewness = trapz(Smooth_AllK, (((R_T - alpha_3_mean) / alpha_3_std_dev) .^ 3) .* P_Density_alpha_3);
alpha_3_kurtosis = trapz(Smooth_AllK, (((R_T - alpha_3_mean) / alpha_3_std_dev) .^ 4) .* P_Density_alpha_3) - 3;


% alpha = 0.04
alpha_4_mean     = trapz(Smooth_AllK, R_T .* P_Density_alpha_4);
alpha_4_variance = trapz(Smooth_AllK, ((R_T - alpha_4_mean) .^ 2) .* P_Density_alpha_4);
alpha_4_std_dev  = sqrt(alpha_4_variance);
alpha_4_skewness = trapz(Smooth_AllK, (((R_T - alpha_4_mean) / alpha_4_std_dev) .^ 3) .* P_Density_alpha_4);
alpha_4_kurtosis = trapz(Smooth_AllK, (((R_T - alpha_4_mean) / alpha_4_std_dev) .^ 4) .* P_Density_alpha_4) - 3;


%%  [Power Utility]

% gamma = 1
gamma_1_mean     = trapz(Smooth_AllK, R_T .* P_Density_gamma_1);
gamma_1_variance = trapz(Smooth_AllK, ((R_T - gamma_1_mean) .^ 2) .* P_Density_gamma_1);
gamma_1_std_dev  = sqrt(gamma_1_variance);
gamma_1_skewness = trapz(Smooth_AllK, (((R_T - gamma_1_mean) / gamma_1_std_dev) .^ 3) .* P_Density_gamma_1);
gamma_1_kurtosis = trapz(Smooth_AllK, (((R_T - gamma_1_mean) / gamma_1_std_dev) .^ 4) .* P_Density_gamma_1) - 3;


% gamma = 2
gamma_2_mean     = trapz(Smooth_AllK, R_T .* P_Density_gamma_2);
gamma_2_variance = trapz(Smooth_AllK, ((R_T - gamma_2_mean) .^ 2) .* P_Density_gamma_2);
gamma_2_std_dev  = sqrt(gamma_2_variance);
gamma_2_skewness = trapz(Smooth_AllK, (((R_T - gamma_2_mean) / gamma_2_std_dev) .^ 3) .* P_Density_gamma_2);
gamma_2_kurtosis = trapz(Smooth_AllK, (((R_T - gamma_2_mean) / gamma_2_std_dev) .^ 4) .* P_Density_gamma_2) - 3;


% gamma = 3
gamma_3_mean     = trapz(Smooth_AllK, R_T .* P_Density_gamma_3);
gamma_3_variance = trapz(Smooth_AllK, ((R_T - gamma_3_mean) .^ 2) .* P_Density_gamma_3);
gamma_3_std_dev  = sqrt(gamma_3_variance);
gamma_3_skewness = trapz(Smooth_AllK, (((R_T - gamma_3_mean) / gamma_3_std_dev) .^ 3) .* P_Density_gamma_3);
gamma_3_kurtosis = trapz(Smooth_AllK, (((R_T - gamma_3_mean) / gamma_3_std_dev) .^ 4) .* P_Density_gamma_3) - 3;


% gamma = 4
gamma_4_mean     = trapz(Smooth_AllK, R_T .* P_Density_gamma_4);
gamma_4_variance = trapz(Smooth_AllK, ((R_T - gamma_4_mean) .^ 2) .* P_Density_gamma_4);
gamma_4_std_dev  = sqrt(gamma_4_variance);
gamma_4_skewness = trapz(Smooth_AllK, (((R_T - gamma_4_mean) / gamma_4_std_dev) .^ 3) .* P_Density_gamma_4);
gamma_4_kurtosis = trapz(Smooth_AllK, (((R_T - gamma_4_mean) / gamma_4_std_dev) .^ 4) .* P_Density_gamma_4) - 3;


%%  Summary

Summary_table = table(...
    [RND_mean; alpha_1_mean; alpha_2_mean; alpha_3_mean; alpha_4_mean; gamma_1_mean; gamma_2_mean; gamma_3_mean; gamma_4_mean], ...
    [RND_variance; alpha_1_variance; alpha_2_variance; alpha_3_variance; alpha_4_variance; gamma_1_variance; gamma_2_variance; gamma_3_variance; gamma_4_variance], ...
    [RND_skewness; alpha_1_skewness; alpha_2_skewness; alpha_3_skewness; alpha_4_skewness; gamma_1_skewness; gamma_2_skewness; gamma_3_skewness; gamma_4_skewness], ...
    [RND_kurtosis; alpha_1_kurtosis; alpha_2_kurtosis; alpha_3_kurtosis; alpha_4_kurtosis; gamma_1_kurtosis; gamma_2_kurtosis; gamma_3_kurtosis; gamma_4_kurtosis], ...
    'VariableNames', {'Mean', 'Variance', 'Skewness', 'Kurtosis'}, ...
    'RowNames', {'RND', 'alpha = 0.01', 'alpha = 0.02', 'alpha = 0.03', 'alpha = 0.04', 'gamma = 1', 'gamma = 2', 'gamma = 3', 'gamma = 4'});

disp(Summary_table);

writetable(Summary_table, fullfile(currentFolder, 'Table', 'Summary_table.csv'), 'WriteRowNames', true);


%%    Summary

Summary_table_T = table(...
    [RND_mean;     RND_variance;     RND_skewness;     RND_kurtosis    ], ...
    [alpha_1_mean; alpha_1_variance; alpha_1_skewness; alpha_1_kurtosis], ...
    [alpha_2_mean; alpha_2_variance; alpha_2_skewness; alpha_2_kurtosis], ...
    [alpha_3_mean; alpha_3_variance; alpha_3_skewness; alpha_3_kurtosis], ...
    [alpha_4_mean; alpha_4_variance; alpha_4_skewness; alpha_4_kurtosis], ...
    [gamma_1_mean; gamma_1_variance; gamma_1_skewness; gamma_1_kurtosis], ...
    [gamma_2_mean; gamma_2_variance; gamma_2_skewness; gamma_2_kurtosis], ...
    [gamma_3_mean; gamma_3_variance; gamma_3_skewness; gamma_3_kurtosis], ...
    [gamma_4_mean; gamma_4_variance; gamma_4_skewness; gamma_4_kurtosis], ...
    'RowNames', {'Mean', 'Variance', 'Skewness', 'Kurtosis'}, ...
    'VariableNames', {'RND', ...
                      'alpha = 0.01', 'alpha = 0.02', 'alpha = 0.03', 'alpha = 0.04', ...
                      'gamma = 1',    'gamma = 2',    'gamma = 3',    'gamma = 4'});

disp(Summary_table_T);

writetable(Summary_table_T, fullfile(currentFolder, 'Table', 'Summary_table_T.csv'), 'WriteRowNames', true);


%% Check

Div = Data(1, Index_DY);
theoretical_value = exp((RF - Div) * (Target_TTM / 365));
