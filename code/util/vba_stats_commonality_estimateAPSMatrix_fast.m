function  APSMatrix = vba_stats_commonality_estimateAPSMatrix_fast(dataTemp,ModelPred,numcc,Ns,y,mu_y)

APSMatrix = nan(numcc, 1);
ss_tot    = sum((y - mu_y).^2);

if ss_tot < eps
    warning('estimateAPSMatrix_fast: response variable has near-zero variance.');
    return
end

% Convert table once and build name?index map (O(1) lookup in loop)
colNames = dataTemp.Properties.VariableNames;
colMap   = containers.Map(colNames, 1:numel(colNames));
dataMat  = table2array(dataTemp);          % [Ns x numPredictors]
intercept = ones(Ns, 1);

for i = 1:numcc

    % Normalise predictor names to cell array
    predNames = cellstr(ModelPred{i});

    % O(1) column lookup via Map ? no ismember scan each iteration
    colIdx = cell2mat(values(colMap, predNames));

    % Calculate R2 by hand
    X      = [intercept, dataMat(:, colIdx)];
    betas  = X \ y;
    ss_res = sum((X * betas - y).^2);
    APSMatrix(i) = 1 - ss_res / ss_tot;
end

% APSMatrix = nan(numcc,1);
% for i = 1:numcc
%     % Calculate R2 by hand
%     X = [ones(Ns,1) dataTemp{:,string(ModelPred{i})}];
%     betas   = X\y;
%     ss_tot  = sum((y - mu_y).^2);
%     ss_res  = sum((X*betas - y).^2);
%     APSMatrix(i) = 1 - ss_res/ss_tot;
% end
