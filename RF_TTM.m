function RF_Target = RF_TTM(Data, Date_Target, TTM_Target)
% Data: [1. Date (YYYYMMDD) | 2. TTM (Days) | 3. RF (Annualized)]

% Number of Unique Pairs
[AllPair, ~, Index_Pair] = unique([Date_Target TTM_Target], 'rows');

for d = 1:size(AllPair, 1)
    % Data on Specific Target Date
    for i = 0:10
        Index = find(Data(:, 1)==str2num(datestr(datenum(num2str(AllPair(d, 1)), 'yyyymmdd') - i, 'yyyymmdd')));      
        if length(Index) > 0
            Data_TTM = Data(Index, 2);
            Data_RF = Data(Index, 3);
        
            break
        else
        end
    end
    clear Index       

    % Interpolate (Linear)
    RF_Pair(d, 1) = interp1(Data_TTM, Data_RF, AllPair(d, 2), 'linear');   
    clear Data_TTM Data_RF
    
%     % Interpolate (Cubic Spline)
%     RF_Pair(d, 1) = interp1(Data_TTM, Data_RF, AllPair(d, 2), 'spline');   
%     clear Data_TTM Data_RF
end
clear AllPair

% Summary
RF_Target = RF_Pair(Index_Pair);
clear Index_Pair RF_Pair
clear d i

end