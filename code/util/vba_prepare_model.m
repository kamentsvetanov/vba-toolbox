function varout = vba_prepare_model(cfg, T)
%
% VBA_PREPARE_MODEL  Parse Wilkinson formula and prepare analysis table.
%
% Handles standard GLM formulas and extended mixed effects formulas.
% For mixed models, parses both fixed and random effects components,
% validates grouping variables, and prepares the full model specification
% so vba_run and estimation subfunctions receive a clean, unambiguous
% model description.
%
% It also, extends model parsing to inspect variable types in tbl and
% dummy-code any categorical variables once, returning a fully
% numeric table ready for QR estimation and commonality analysis.
%
% INPUTS
%   Model - Wilkinson model string. Two formats accepted:
%
%           GLM (linear/commonality):
%             'f_rsfa ~ age + sex'
%
%           Mixed effects ? random effects appended in parentheses:
%             'f_rsfa ~ age + sex + (1|SubjectID)'
%             'f_rsfa ~ age + time + (1 + time|SubjectID)'
%             'f_rsfa ~ age + time + (1|SubjectID) + (1|Site)'
%
%   tbl   - Subject table (optional). If provided:
%             - categorical variables are dummy-coded
%             - grouping variables are validated
%             - random effect levels are counted and reported
%
%   cfg   - Config struct (optional). Uses cfg.modelType to determine
%           whether to expect random effects syntax.
%
% OUTPUTS
%   varout.ResponseVar       - Cell, response variable name
%   varout.PredictorVars     - Cell, fixed effect predictor names (no random)
%   varout.VarNames          - Cell, all variable names (expanded for dummies)
%   varout.f_model           - String, sanitised model name for filenames
%   varout.moderation        - Struct, interaction terms if present
%   varout.fixedFormula      - String, fixed effects only formula for fitlm
%   varout.mixedFormula      - String, full formula for fitlme (if mixed)
%   varout.randomEffects     - Struct array, one per random effects term:
%                                .formula    '(1 + time|SubjectID)'
%                                .slopes     {'time'} or {} for intercept only
%                                .groupVar   'SubjectID'
%                                .numGroups  number of unique levels (if tbl given)
%                                .intercept  logical, random intercept included
%   varout.tbl               - Table, categoricals dummy-coded (if tbl given)
%   varout.colMap            - containers.Map, varname->column indices
%   varout.dummyInfo         - Struct array, dummy coding reference
%
% Author : Kamen Tsvetanov, Ph.D.
%__________________________________________________________________________

modelType = '';
if ~isempty(cfg) && isfield(cfg,'modelType')
    modelType = cfg.modelType;
end

Model = cfg.model;

% =========================================================================
% STEP 1 - Split fixed and random effects  (unchanged)
% =========================================================================
[fixedFormula, randomEffects] = parseRandomEffects(Model);

% =========================================================================
% STEP 2 - Parse fixed effects formula  (unchanged)
% =========================================================================
VarNames      = strtrim(split(fixedFormula, ["~","*","+","^",":"]));
VarNames      = VarNames(~cellfun(@isempty, VarNames));
ResponseVar   = VarNames(1);
PredictorVars = VarNames(2:end);

% =========================================================================
% STEP 3 - Handle intercept-only / main effect model
% '1' as sole predictor means one-sample t-test; suppress z-scoring
% =========================================================================
isMainEffect = any(strcmp(PredictorVars, '1'));
doZscore     = cfg.doZscore;   % may be overridden below

if isMainEffect
    PredictorVars = PredictorVars(~strcmp(PredictorVars, '1'));
    VarNames      = VarNames(~strcmp(VarNames, '1'));
    doZscore      = 0;
    fprintf('  Main effect model detected: z-scoring disabled\n');
end

% =========================================================================
% STEP 4 - Resolve covariate prefix 'c_' in T
%
% Variables named 'c_xxx' in the model mark covariates of no interest
% (excluded from output maps). If T does not contain 'c_xxx' as a column
% but does contain 'xxx', create 'c_xxx' as an alias so downstream
% code can find it by its prefixed name.
% =========================================================================
idxCovariate  = contains(VarNames, 'c_');
nameCovariate = VarNames(idxCovariate);
VarNamesTable = T.Properties.VariableNames;

for ivar = 1:numel(nameCovariate)
    namecov         = nameCovariate{ivar};
    namecovNoPrefix = namecov(3:end);                      % strip 'c_'

    hasPrefix   = any(ismember(VarNamesTable, namecov));
    hasNoPrefix = any(ismember(VarNamesTable, namecovNoPrefix));

    if ~hasPrefix && hasNoPrefix
        % Alias: add prefixed column pointing to existing data
        T.(namecov) = T.(namecovNoPrefix);
        fprintf('  Covariate alias: ''%s'' -> ''%s''\n', namecov, namecovNoPrefix);
    elseif ~hasPrefix && ~hasNoPrefix
        error(['vba_prepare_model: covariate ''%s'' (or ''%s'') not found in T.\n' ...
               'Available variables: %s'], ...
               namecov, namecovNoPrefix, strjoin(VarNamesTable, ', '));
    end
end

% =========================================================================
% STEP 5 - Filter subjects with missing data
%
% Only check missingness for non-imaging variables here ? imaging
% variables (f_ prefix) are checked per-voxel in vba_load_images
% via the valid-voxel mask. Including f_ columns here would require
% loading all images just to filter subjects, which is premature.
% =========================================================================
VarNamesNonImaging = VarNames(~contains(VarNames, 'f_'));

if ~isempty(VarNamesNonImaging)
    % Keep only columns that exist in T (imaging cols may not yet)
    VarNamesCheck = VarNamesNonImaging(ismember(VarNamesNonImaging, VarNamesTable));
    idxSub        = all(~ismissing(T(:, VarNamesCheck)), 2);
else
    idxSub        = true(height(T), 1);
end

numExcluded = sum(~idxSub);
if numExcluded > 0
    fprintf('  Excluded %d subject(s) with missing non-imaging data\n', numExcluded);
end

tbl = T(idxSub, :);   % full row subset, all columns retained

% =========================================================================
% STEP 6 - Validate random effects against tbl  (unchanged)
% =========================================================================
if ~isempty(randomEffects) && ~isempty(tbl)
    randomEffects = validateRandomEffects(randomEffects, tbl);
end

% =========================================================================
% STEP 7 - Dummy code categoricals in tbl  (unchanged)
% =========================================================================

[tbl, colMap, dummyInfo, PredictorVars, VarNames] = ...
    dummyCodeTable(tbl, PredictorVars, VarNames);

% =========================================================================
% STEP 8 - Reconstruct formulas  (unchanged)
% =========================================================================
% fixedFormula  - for fitlm or the QR path (no random effects)
% mixedFormula  - for fitlme (fixed + random)
fixedFormulaClean = rebuildFixedFormula(ResponseVar{1}, PredictorVars, isMainEffect);

if ~isempty(randomEffects)
    reParts      = {randomEffects.formula};
    mixedFormula = [fixedFormulaClean ' + ' strjoin(reParts, ' + ')];
else
    mixedFormula = fixedFormulaClean;
end

% =========================================================================
% Assemble output
% =========================================================================
varout.ResponseVar      = ResponseVar;
varout.PredictorVars    = PredictorVars;
varout.VarNames         = VarNames;
varout.f_model          = sanitiseModelName(fixedFormula);
varout.isMainEffect     = isMainEffect;
varout.doZscore         = doZscore;        % may differ from cfg.doZscore
varout.fixedFormula     = fixedFormulaClean;
varout.mixedFormula     = mixedFormula;
varout.randomEffects    = randomEffects;
varout.tbl              = tbl;
varout.colMap           = colMap;
varout.dummyInfo        = dummyInfo;
varout.idxSub           = idxSub;          % logical index into original T


% Print summary
printSummary(varout, cfg.modelType);

end


% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function [fixedFormula, randomEffects] = parseRandomEffects(Model)
%
% Extract random effects terms (slope|group) from formula.
% Returns the fixed-effects-only formula and a struct array
% describing each random effects term.

% Find all parenthetical random effects terms: (... | ...)
rePattern = '\(([^)]+)\)';
reTokens  = regexp(Model, rePattern, 'tokens');

randomEffects = struct('formula',{},'slopes',{},'groupVar',{},...
                       'intercept',{},'numGroups',{});

fixedFormula = Model;

for ir = 1:numel(reTokens)
    inner   = strtrim(reTokens{ir}{1});    % e.g. '1 + time|SubjectID'
    fullTerm= ['(' inner ')'];

    % Split on '|' to get slopes | grouping variable
    parts = strsplit(inner, '|');
    if numel(parts) ~= 2
        warning('vba_prepare_model: cannot parse random effects term ''%s'', skipping.', fullTerm);
        continue
    end

    slopePart = strtrim(parts{1});         % e.g. '1 + time' or '1'
    groupVar  = strtrim(parts{2});         % e.g. 'SubjectID'

    % Parse slope components
    slopeTerms   = strtrim(split(slopePart, '+'));
    hasIntercept = any(strcmp(slopeTerms, '1'));
    slopes       = slopeTerms(~strcmp(slopeTerms,'1'));  % non-intercept slopes

    % Populate struct
    re.formula   = fullTerm;
    re.slopes    = slopes;
    re.groupVar  = groupVar;
    re.intercept = hasIntercept;
    re.numGroups = [];          % filled later if tbl is provided
    randomEffects(end+1) = re; %#ok<AGROW>

    % Remove this term from the fixed formula
    % Also remove any leading/trailing ' + ' left behind
    fixedFormula = strrep(fixedFormula, [' + ' fullTerm], '');
    fixedFormula = strrep(fixedFormula, [fullTerm ' + '], '');
    fixedFormula = strrep(fixedFormula, fullTerm, '');
    fixedFormula = strtrim(fixedFormula);
end
end


function randomEffects = validateRandomEffects(randomEffects, tbl)
%
% Check grouping variables exist in tbl and count their levels.
% Warns about groups with very few levels (< 5) which can cause
% singular fits in mixed models.

colNames = tbl.Properties.VariableNames;

for ir = 1:numel(randomEffects)
    gv = randomEffects(ir).groupVar;

    if ~ismember(gv, colNames)
        error(['vba_prepare_model: grouping variable ''%s'' not found in tbl.\n' ...
               'Available variables: %s'], gv, strjoin(colNames, ', '));
    end

    levels    = unique(tbl.(gv));
    numGroups = numel(levels);
    randomEffects(ir).numGroups = numGroups;

    fprintf('  Random effect (1|%s): %d groups\n', gv, numGroups);

    if numGroups < 5
        warning(['vba_prepare_model: grouping variable ''%s'' has only %d levels.\n' ...
                 'Mixed models typically require more groups for stable variance estimates.'], ...
                 gv, numGroups);
    end

    % Validate random slopes exist as fixed effects too
    for is = 1:numel(randomEffects(ir).slopes)
        slope = randomEffects(ir).slopes{is};
        if ~ismember(slope, tbl.Properties.VariableNames)
            error(['vba_prepare_model: random slope ''%s'' in term ''%s'' ' ...
                   'not found in tbl.'], slope, randomEffects(ir).formula);
        end
    end
end
end


function [tbl, colMap, dummyInfo, PredictorVars, VarNames] = ...
         dummyCodeTable(tbl, PredictorVars, VarNames)
%
% Inspect non-imaging predictors in tbl.
% Numeric columns pass through; categoricals are dummy-coded
% with the first level as reference (matching fitlm behaviour).

colMap    = containers.Map();
dummyInfo = struct('varName',{},'levels',{},'refLevel',{},'dummyNames',{});

nonImagingPreds = PredictorVars(~contains(PredictorVars,'f_'));

for iv = 1:numel(nonImagingPreds)
    varName = nonImagingPreds{iv};

    if ~ismember(varName, tbl.Properties.VariableNames)
        warning('vba_prepare_model: ''%s'' not found in tbl, skipping.', varName);
        continue
    end

    v = tbl.(varName);

    if isnumeric(v) || islogical(v)
        colMap(varName) = varName;  % resolves to column name; index set later

    elseif iscategorical(v) || iscell(v) || ischar(v) || isstring(v)
        [idx, levels] = grp2idx(v);
        numLevels     = numel(levels);

        if numLevels <= 1
            warning('vba_prepare_model: ''%s'' has only one level, skipping.', varName);
            continue
        end

        % Create dummy columns (drop first level as reference)
        dummyNames = cell(1, numLevels-1);
        for lv = 2:numLevels
            dname           = sprintf('%s_%s', varName, string(levels(lv)));
            dummyNames{lv-1}= dname;
            tbl.(dname)     = double(idx == lv);
        end

        % Remove original column
        tbl.(varName) = [];

        % Log
        di                   = numel(dummyInfo) + 1;
        dummyInfo(di).varName    = varName;
        dummyInfo(di).levels     = levels;
        dummyInfo(di).refLevel   = levels(1);
        dummyInfo(di).dummyNames = dummyNames;

        % Map original name to its dummy names
        colMap(varName) = dummyNames;

        % Expand PredictorVars and VarNames
        PredictorVars = [PredictorVars(~strcmp(PredictorVars,varName)), dummyNames];
        % VarNames      = [VarNames(~strcmp(VarNames,varName)),           dummyNames];
        VarNames = [reshape(VarNames(~strcmp(VarNames,varName)), 1, []), ...
            reshape(dummyNames, 1, [])];

        fprintf('  Dummy-coded ''%s'': ref=''%s'', %d levels -> %d columns\n', ...
            varName, string(levels(1)), numLevels, numLevels-1);
    else
        warning('vba_prepare_model: ''%s'' has unsupported type, skipping.', varName);
    end
end

% Resolve all colMap entries to numeric column indices
allCols = tbl.Properties.VariableNames;
keys_   = keys(colMap);
for ik = 1:numel(keys_)
    k = keys_{ik};
    v = colMap(k);
    if ischar(v) || isstring(v)
        colMap(k) = find(strcmp(allCols, v));
    elseif iscell(v)
        colMap(k) = cellfun(@(d) find(strcmp(allCols,d)), v);
    end
end
end


function formula = rebuildFixedFormula(responseVar, PredictorVars, isMainEffect)
%
% Reconstruct a clean fixed-effects Wilkinson formula string
% from the parsed (and possibly dummy-expanded) components.

if isMainEffect || isempty(PredictorVars)
    formula = sprintf('%s ~ 1', responseVar);
else
    formula = sprintf('%s ~ %s', responseVar, strjoin(PredictorVars, ' + '));
end
end


function moderation = parseInteractions(formula, ResponseVar, PredictorVars)
%
% Extract interaction terms from fixed effects formula.

moderation.Y = ResponseVar{1};
str = regexp(formula,'(\w*)[*|:](\w*)[*|:](\w*)','tokens');
if ~isempty(str)
    moderation.X = str{1}{1};
    moderation.W = str{1}{2};
    moderation.Z = str{1}{3};
    moderation.C = PredictorVars(~ismember(PredictorVars, ...
                   [moderation.X, moderation.W, moderation.Z]));
else
    str = regexp(formula,'(\w*)[*|:](\w*)','tokens');
    if ~isempty(str)
        moderation.X = str{1}{1};
        moderation.Z = str{1}{2};
        moderation.C = PredictorVars(~ismember(PredictorVars, ...
                       [moderation.X, moderation.Z]));
    end
end
end


function printSummary(varout, modelType)
fprintf('\n--- Model specification ---\n');
fprintf('  Fixed formula : %s\n', varout.fixedFormula);
if ~isempty(varout.randomEffects)
    fprintf('  Mixed formula : %s\n', varout.mixedFormula);
    for ir = 1:numel(varout.randomEffects)
        re = varout.randomEffects(ir);
        if isempty(re.slopes)
            slopeStr = 'intercept only';
        else
            slopeStr = ['slopes: ' strjoin(re.slopes,', ')];
        end
        fprintf('    RE term %d   : %s | group=%s (%s)', ...
            ir, re.formula, re.groupVar, slopeStr);
        if ~isempty(re.numGroups)
            fprintf(', %d groups', re.numGroups);
        end
        fprintf('\n');
    end
end
fprintf('  Response      : %s\n', varout.ResponseVar{1});
fprintf('  Predictors    : %s\n', strjoin(varout.PredictorVars, ', '));
if ~isempty(modelType)
    fprintf('  Model type    : %s\n', modelType);
end
fprintf('---------------------------\n\n');
end


function strModel = sanitiseModelName(fixedFormula)
    strModel = fixedFormula;
    strModel = regexprep(strModel, 'f_',  '');
    strModel = regexprep(strModel, 'c_',  '');
    strModel = regexprep(strModel, '(\<[a-z])','${upper($1)}');
    strModel = regexprep(strModel, {'+',' '},'');
    strModel = regexprep(strModel, '~',  '_iv');
    strModel = regexprep(strModel, '*',  'X');
    strModel = ['dv' strModel];
end