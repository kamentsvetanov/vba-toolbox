function tvals = vba_filtm(cfg, datY, iperm, tbl, VarNamesMaps, parforArg)

randOrder   = cfg.randOrder;
doRobust    = cfg.doRobust;
numVox      = cfg.numVox;
numCoef     = cfg.numCoef;
ResponseVar = cfg.mlr.ResponseName;
model       = cfg.model;

tvals       = zeros(numVox,numCoef);
tempOrder   = randOrder(:,iperm);

% Loop through all voxels
 parfor (iVox = 1:numVox, parforArg) 
   
    Yvox = squeeze(datY(:,iVox,:));

    % Ensure all subjects have non-nan values
    if sum(sum(isnan(Yvox)))==0
        tbl_temp = tbl;    
        tbl_temp{:,VarNamesMaps} = Yvox;
        tbl_temp(:,ResponseVar) = tbl_temp(tempOrder,ResponseVar);% Shuffle dependent variable

        mlr             = [];  
        mlr             = fitlm(tbl_temp,model,'RobustOpts',doRobust);
%                 bvals(iVox,:)   = mlr.Coefficients.Estimate;
        tvals(iVox,:)   = mlr.Coefficients.tStat;
    end
end

