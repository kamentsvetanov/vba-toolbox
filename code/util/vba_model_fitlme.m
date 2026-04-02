function tvals = vba_model_fitlme(cfg, Y, iperm, tbl, VarNamesMaps, parforArg)

numVox  = cfg.numVox;
numCoef = cfg.numCoef;
Model   = cfg.model;
randOrder = cfg.randOrder;
tvals   = nan(numCoef, numVox);


randOrder   = cfg.randOrder;
doRobust    = cfg.doRobust;
numVox      = cfg.numVox;
numCoef     = cfg.numCoef;
ResponseVar = cfg.mlr.ResponseName;
model       = cfg.model;

tvals       = zeros(numVox,numCoef);
tempOrder   = randOrder(:,iperm);
% Permute the DV rows for this permutation
permIdx = randOrder(:, iperm);

parfor (ivox = 1:numVox, parforArg)
    tbl_vox = tbl;

    % Inject imaging data for this voxel
    for iMap = 1:numel(VarNamesMaps)
        tbl_vox.(VarNamesMaps{iMap}) = Y(permIdx, ivox, iMap);
    end

    % Skip voxels with insufficient variance or all-NaN
    dvName = cfg.ResponseVar;
    if var(tbl_vox.(dvName), 'omitnan') < eps
        continue
    end

    try
        lme = fitlme(tbl_vox, Model, ...
            'FitMethod',    'REML', ...   % use ML if comparing fixed effects
            'Verbose',      false, ...
            'CheckHessian', false);       % skip for speed; re-enable for diagnostics

        if lme.converged
            tvals(:, ivox) = lme.Coefficients.tStat;
        end
        % Could also save: lme.Coefficients.pValue, lme.Coefficients.Estimate

    catch
        % Convergence failure or rank deficiency ? leave as NaN
    end
end
end