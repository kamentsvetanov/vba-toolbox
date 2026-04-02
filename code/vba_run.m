function [cfg,outDir] = vba_run(T,cfg) 
%
% A function for group-level multi-modal voxel-wise brain image analysis.
% Running the analysis requires two variables in the workspace, T and
% Model, similar to the usage of 'fitlm'
% 
% T       - Table containing all variables specified in Model. Similar to
%           'fitlm', variables in T can be continuous or categorical
%           variables. For voxel-based analysis is required to specify
%           filepaths to individuals' images (e.g. anatomical T1w scans) in
%           a variable prefixed by 'f_', e.g. f_t1w for a structure
%           variable in T containing filepaths to anatomical scans.
%
% Model   - A string specifying the linear model formula using Wilkinson notation.
%           (https://uk.mathworks.com/help/stats/wilkinson-notation.html)
%           y ~ terms, where y is the name of the response variable, and
%           terms defines the model using the predictor variable names and
%           the operators.
%           Images (e.g. anatomical T1w scans) can be used as response and/
%           or predictor variables. 
%           Model = 'f_t1w ~ age + sex'; Specifies analysis to predict T1w
%           intensity using age and sex as predictors.
% 
% This GLM-like approach uses fitlm to save nii image of coeff and p-values 
% for each variable and residuals (unexplained effects by predictors).
% This could be useful in instances with voxel-specific covariates, e.g. to
% estimate variance explained or residuals in RSFA maps (across subjects)
% after controlling for the effects of ASL maps (voxel-specific), HRV and
% other effects as in 
% Tsvetanov et al 2020 Psyhophysiology (https://doi.org/10.1111/psyp.13714)
% and
% Wu et al 2021 (https://www.biorxiv.org/content/10.1101/2021.11.10.468042v1)
% 
% Notes:
% Commonality Analysis:
% - voxel-based analysis (involving thousands of voxels) is slow useing
%   'fitlm'. A more efficient version is implemented when non-robust glm is
%   acceptable. Note that this approach, unlike fitlm, works only with
%   continuous variables i.e. any categorical variable needs to be used as
%   continous (e.g. use grp2idx). Categorical variables can be excluded
%   from the output by naming it with 'c_' prefix. 
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2
%
% Script edit log (kat)
% ---------------------
% 02/10/2017 - function initiation
% 12/10/2021 - streamlining with WIlkinson notations
% 10/01/2022 - estimate main effects as 'One-Sample' by specifying 'y ~ 1'
% 15/04/2022 - choose optional starting value of the permutation number
%              (e.g. continue aborted analysis)
% 21/04/2022 - omit reporting covariates of no interest by prefixing with 'c_'      
% 30/05/2024 - PALM library (palm_quickperms) no longer required
%
% Author : Kamen Tsvetanov, Ph.D., Neurocognitive Ageing
% Affil. : Department of Clinical Neurosciences, University of Cambridge
% Email  : kamen.tsvetanov@gmail.com  
% Website: http://www.kamentsvetanov.com
%__________________________________________________________________________
% Copyright (C) Kamen Tsvetanov 2018


% try rootDir         = cfg.rootDir;          catch rootDir       = pwd;  end % Root directory to output analysis results
% try doZscore        = cfg.doZscore;         catch doZscore      = 1;    end % Whether or not to z-score data (Default,1)
% try doRobust        = cfg.doRobust;         catch doRobust      = 0;    end % Whether or not to run robust regression (Default, 0)
% try doRunInSerial   = cfg.doRunInSerial;    catch doRunInSerial = 1;    end % whether or not to run in Parallel mode
% try numPerm         = cfg.numPerm;          catch numPerm       = 200;  end % Number of permuations
% try startPerm       = cfg.startPerm;        catch startPerm     = 1;    end % (Optional) pick a different starting position of the permuations
% try f_mask          = cfg.f_mask;           catch error('Error. \nPlease provide brain mask.'); end % Brain mask to limit analysis for voxels in the mask
% try Model           = cfg.model;            catch error('Error. \nPlease specify the model using Wilkinson notaions.'); end 
% try doSlurm         = cfg.doSlurm;          catch doSlurm         = 0;  end % (SLURM option)
% try predefRandOrder = cfg.predefRandOrder;  catch predefRandOrder = []; end % (SLURM option) provde a predefined permuted order of DV 
% try whichRandOrder  = cfg.whichRandOrder;   catch whichRandOrder  = []; end % (SLURM option) specify which order to select from the permuted matrix
% try specificSeed    = cfg.specificSeed;     catch specificSeed    = []; end % (SLURM option) provide specificSeed to reproduce permuted Matrix across workers/nodes
% try doLogTrans      = cfg.doLogTrans;       catch doLogTrans    = 0;    end % Whether or not to log-transform neuroimaging data

% -------------------------------------------------------------------------
% 1. Validate and fill defaults ? single source of truth
cfg = vba_cfg(cfg);   % populate defaults and validate

% -------------------------------------------------------------------------
% 2. Unpack for readability and parfor safety
% Scalars and strings only ? large arrays stay in cfg
% and are passed explicitly to subfunctions that need them
rootDir          = cfg.rootDir;
doZscore         = cfg.doZscore;
doLogTrans       = cfg.doLogTrans;
doRobust         = cfg.doRobust;
doRunInSerial    = cfg.doRunInSerial;
numPerm          = cfg.numPerm;
startPerm        = cfg.startPerm;
Model            = cfg.model;
modelType        = cfg.modelType;% This can be 'linear', 'mixed' or 'commonality' or 'maineffect'
specificSeed     = cfg.specificSeed;
whichRandOrder   = cfg.whichRandOrder;
f_mask           = cfg.f_mask;



% -------------------------------------------------------------------------
% Set parforArg to inf (Default), assuming CA is run in matlab with parfor.
% Otherwise set to 0 (see example with SLURM implementation below)
parforArg = Inf;
if doRunInSerial % Flag Paraller processing
    parforArg = 0;
end


% modeltype       = 'regression'; % 'one-sample' also possible if 'y~1'
idxNumeric      = vartype('numeric'); 
% idxCategorical  = vartype('categorical'); 

% -------------------------
% Get Variable names
% -------------------------
M               = vba_prepare_model(cfg,T);
% Unpack everything downstream needs
tbl           = M.tbl;
VarNames      = M.VarNames;
PredictorVars = M.PredictorVars;
ResponseVar   = M.ResponseVar;
fixedFormula  = M.fixedFormula;
mixedFormula  = M.mixedFormula;   % passed to vba_model_fitlme
randomEffects = M.randomEffects;  % passed into cfg for mixed subfunction
cfg.colMap    = M.colMap;
cfg.doZscore  = M.doZscore;     % may have been forced to 0 by main effect
Model         = M.fixedFormula;
cfg.model     = Model;


%----------------------------
% Load imaging data 
% ---------------------------
[Y, meta] = vba_load_images(tbl, VarNames, cfg);

% Unpack what vba_run needs directly
Vdv             = meta.Vdv;
idxMask         = meta.idxMask;
idxMaskValid    = meta.idxMaskValid;
VarNamesMaps    = meta.VarNamesMaps;
VarNamesGlobal  = meta.VarNamesGlobal;
numVox          = meta.numVox;
numSub          = meta.numSub;

% Update cfg with resolved sizing
cfg.numVox      = meta.numVox;
cfg.numSub      = meta.numSub;
cfg.idxValidVox = meta.idxValidVox;



% %----------------------------
% % Import and apply Brain mask 
% % ---------------------------
% Ymask   = logical(spm_read_vols(spm_vol(f_mask)));
% dimMask = size(Ymask);
% idxMask = find(Ymask);
% 
% % --------------------------------------------------------------
% % Import Variables/Modalities with Brain images (e.g. RSFA maps
% % ---------------------------------------------------------------
% idxMaps         = find(contains(VarNames,'f_'));
% VarNamesMaps    = VarNames(idxMaps);
% VarNamesGlobal  = VarNames(~contains(VarNames,'f_'));
% 
% numVox   = size(idxMask,1);
% numSub   = size(tbl,1);
% numMap   = numel(idxMaps);
% 
% Y = nan(numSub,numVox,numMap);
% 
% for iMap = 1:numMap
%     nameMaps = VarNames{idxMaps(iMap)};
%     V   = spm_vol(char(tbl.(nameMaps)));
%     % Check that Mask and this modality dimensions match
%     if ~isequal(dimMask,V(1).dim)
%         strError = sprintf('Error. \nDimensions of Mask and %s do not match.',nameMaps);
%         error(strError);
%     end       
%     y   = spm_read_vols(V);
%     y   = permute(y,[4 1 2 3]);
%     Y(:,:,iMap) = y(:,Ymask);
% 
%     if iMap == 1
%         Vdv = V(1); % Store head of dependent variable modality for later
%     end
% end
% 
% % Log-transform imaging data
% % --------------------------
% % N.B. needs further refinement, as it would apply to all modalities (if
% % more than one modality is given).
% if doLogTrans
%     Y = log(Y);
% end

tbl_to_save  = tbl;

if doZscore
    tbl(:,idxNumeric) = normalize(tbl(:,idxNumeric));
    Y = normalize(Y,1);
end

% ------------------------------------------------------------------------
% Prepare the table for fitlm with one voxel for dependent and predictor
% maps and set a template to extract data after parpool
% ------------------------------------------------------------------------

Temp = tbl(:,VarNamesGlobal);
Temp{:,VarNamesMaps} = squeeze(Y(:,ceil(numVox/2),:)); %If may crash if the selected voxel has NaN or other issues.Could replace with random variable having same lenght
tbl = Temp;
mlr_temp     = fitlm(tbl,Model);
nameCoef     = mlr_temp.CoefficientNames;


% Clean up workspace to imporve parpool performance
clear T X X3D Xvec Y3D Temp Yvec; 

%% make output folders for each  permutation
strModel = M.f_model;
strModel = [strModel '_n' num2str(numSub) '_nPerm' num2str(numPerm)];
if doRobust
    strModel = [strModel '_Robust' num2str(doRobust)];
end
strModel = [strModel '_' datestr(now,'yyyymmdd')];

% ----------------------------------------------
% Predefine Commonality Analysis or GLM output maps
% ----------------------------------------------
switch modelType
    case 'commonality'
        outDir      = fullfile(rootDir,['ca_' strModel]);
        cfg.mlr     = mlr_temp;
        cfg.doPerm  = 0;
        CA_temp     = vba_stats_commonality(cfg);
        nameVarCA   = CA_temp.Properties.RowNames;
        nameVarCA   = regexprep(nameVarCA,',','');
        numVarCA    = numel(nameVarCA);    
        nameCoef    = nameVarCA;
        numCoef     = numVarCA;
    case 'linear' 
        outDir = fullfile(rootDir,['lm_' strModel]);
        numCoef      = numel(nameCoef);
    case 'mixed' 
        outDir = fullfile(rootDir,['lme_' strModel]);
        numCoef      = numel(nameCoef);
    case 'maineffect' 
        outDir = fullfile(rootDir,['maineffect_' strModel]);
        numCoef      = numel(nameCoef);
end



% -----------------
% Make folders
% -----------------
for i=1:numCoef
    
    
    outname = regexprep(nameCoef{i},'f_','');
    % -------------------------------------------------------------------------
    % Rename interaction terms, so that ':' in interaction terms is replaced 
    % with 'X'. This is for ':' is invalid character for variable and filenames
    % -------------------------------------------------------------------------
    outname = regexprep(outname,':','X');  

    % -------------------------------------------------------------------------
    % Rename squared terms, so that '^' is removed, avoiding invalid character for variable and filenames
    % -------------------------------------------------------------------------
    outname = regexprep(outname,'\^','');  
    nameOutput{i} = outname;

    % Do not create folders for for a pattern specified in excludeFromResult
    if ~contains(nameCoef{i},cfg.excludeFromResults)
        mkdir(fullfile(outDir,['tval_',outname]));
    %     mkdir(fullfile(outDir,['bval_',nameOutput]));
    end
end


%% ------------------------------------------------------------------------
% If SLURM implemenation is required, i.e. predefinedSeed and whichOrder 
% provided, then generate set the seed so that the permuted version 
% using arrayfun (Add 10% extra for now) are identical across workers

if ~isempty(specificSeed)
    rng(specificSeed)
end


% Alternative to palm_quickperms
permutedMatrix = arrayfun(@(x) randperm(numSub), 1:numPerm*1.1, 'UniformOutput', false);
permutedMatrix = cell2mat(permutedMatrix')';
% Remove columns that cointain the original order
idx = corr(permutedMatrix,[1:numSub]')==1;
permutedMatrix(:,idx)=[];
permutedMatrix(:,1) = [1:numSub]'; % Set first column to original order
% get the right size of pertmuted Matrix
randOrder = permutedMatrix(:,1:numPerm);
randOrderAll = randOrder;
% randOrder   = palm_quickperms(numSub,[],numPerm);


if ~isempty(whichRandOrder)
    % Specific randOrder provided. Likely for SLURM impmementation.
    % So reset numPerm to 1 and parpool argument to 0
    randOrder   = randOrder(:,whichRandOrder);
    numPerm     = 1;
    parforArg   = 0;
end


%-Assemble and save output structure
% ----------------------------------
cfg.tbl         = tbl_to_save;
cfg.randOrder   = randOrder;
cfg.randOrderAll= randOrderAll;
cfg.outDir      = outDir;
cfg.mlr         = mlr_temp;
cfg.doZscore    = doZscore;
cfg.doRobust    = doRobust;
cfg.startPerm   = startPerm;
cfg.numVox      = numVox;
cfg.numCoef     = numCoef;

fout = fullfile(outDir,'analysis_cfg.mat');

%if isempty(whichRandOrder) || whichRandOrder==1
% if isempty(whichRandOrder)
    % fout = fullfile(outDir,'analysis_cfg.mat');
% else
%     fout = fullfile(outDir,sprintf('analysis_cfg_%05d.mat',whichRandOrder));
% end
if ~exist(fout,"file") % Check if file has already been generated, e.g. as in SLURM execution
    save(fout,'cfg');
end


for iperm = startPerm:numPerm

    % datY        = Y;
    % --------------------
    % Now run the VBA for the model type
    %---------------------
    switch modelType
        case 'maineffect' % INCOMPLETE. Estimate Main/Average/Group effect using one-sample t-test
            cfg_temp = cfg;
            tvals = vba_model_maineffect(cfg_temp);
        case 'commonality'
            tvals = vba_model_commonality(cfg, Y, iperm, tbl, VarNamesMaps, parforArg);
        case 'linear'
            tvals = vba_model_fitlm(cfg, Y, iperm, tbl, VarNamesMaps, parforArg);
        case 'mixed'  % INCOMPLETE.
            tvals = vba_model_fitlme(cfg,Y, iperm, tbl, VarNamesMaps, parforArg);
    end
       

    % --------------------------------------------------------
    % Write results (coeff, pvals and residuals) to nii images
    % --------------------------------------------------------
    % For every permutation Write results (tvalCA) to nii images in separate folder,
    % so for 1000 permutations there should be 1000 folders
    % (perm_00001,perm_00002 etc), each one containing the tvalCA maps for
    % every effect (common and shared). perm_00001 is the original
    % ordering, i.e. real results.
    if isempty(whichRandOrder)
        vba_write_results(cfg, tvals, nameOutput, Vdv, idxMask, iperm);
    else
        vba_write_results(cfg, tvals, nameOutput, Vdv, idxMask, whichRandOrder);
    end
end



