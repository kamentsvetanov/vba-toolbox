function cfg = vba_cfg(cfg)
%
% VBA_CFG  Set default configuration parameters for voxel-based analysis.
%
% Populates missing fields in cfg with sensible defaults. Any field
% already set by the user is left untouched. Required fields (f_mask,
% model, modelType) throw an informative error if absent.
%
% USAGE
%   cfg = vba_cfg(cfg)       % fill missing fields with defaults
%   cfg = vba_cfg(struct())  % start from scratch with all defaults
%   cfg = vba_cfg()          % same, no input needed
%
% REQUIRED FIELDS (no default ? must be set by user)
%   cfg.f_mask      - Filepath to brain mask NIfTI image
%   cfg.model       - Wilkinson model string, e.g. 'f_rsfa ~ age + sex'
%   cfg.modelType   - 'linear' | 'mixed' | 'commonality' | 'maineffect'
%
% OUTPUT
%   cfg - Configuration struct with all fields populated
%
% Author : Kamen Tsvetanov, Ph.D.
% Affil. : Department of Clinical Neurosciences, University of Cambridge
%__________________________________________________________________________

if nargin < 1
    cfg = struct();
end

% =========================================================================
% REQUIRED ? error if absent
% =========================================================================
cfg = require(cfg, 'f_mask',    'Please provide a brain mask filepath: cfg.f_mask');
cfg = require(cfg, 'model',     'Please specify the model: cfg.model (Wilkinson notation)');
cfg = require(cfg, 'modelType', 'Please specify cfg.modelType: ''linear'', ''mixed'', ''commonality'', or ''maineffect''');

% =========================================================================
% PATHS AND OUTPUT
% =========================================================================
cfg = default(cfg, 'rootDir',   pwd);   % Root directory for output

% =========================================================================
% DATA PROCESSING
% =========================================================================
cfg = default(cfg, 'doZscore',  1);     % Z-score continuous variables
cfg = default(cfg, 'doLogTrans',0);     % Log-transform imaging data

% =========================================================================
% MODEL ESTIMATION
% =========================================================================
cfg = default(cfg, 'doRobust',          0);      % Robust regression (fitlm path only)
cfg = default(cfg, 'estimationMethod',  'QR');   % 'QR' (fast) or 'fitlm' (full)

% Mixed effects specific
cfg = default(cfg, 'lme_randomEffects',     '(1|SubjectID)'); % Random effects formula
cfg = default(cfg, 'lme_groupVar',          'SubjectID');     % Grouping variable name
cfg = default(cfg, 'lme_fitMethod',         'REML');          % 'REML' or 'ML' for fitlme

% =========================================================================
% PERMUTATION TESTING
% =========================================================================
cfg = default(cfg, 'numPerm',   200);   % Number of permutations
cfg = default(cfg, 'startPerm', 1);     % Starting permutation (resume support)

% =========================================================================
% PARALLELISATION
% =========================================================================
cfg = default(cfg, 'doRunInSerial', 1); % 0 = parfor, 1 = serial

% =========================================================================
% SLURM / HPC
% =========================================================================
cfg = default(cfg, 'doSlurm',          0);  % SLURM execution mode
cfg = default(cfg, 'specificSeed',     []); % Seed for reproducible permutations
cfg = default(cfg, 'predefRandOrder',  []); % Predefined permutation matrix
cfg = default(cfg, 'whichRandOrder',   []); % Which column of permutation matrix to use

% =========================================================================
% Validate field values where sensible
% =========================================================================
validModelTypes = {'linear','mixed','commonality','maineffect'};
if ~ismember(cfg.modelType, validModelTypes)
    error('vba_cfg: cfg.modelType ''%s'' is not valid. Choose from: %s.', ...
        cfg.modelType, strjoin(validModelTypes, ', '));
end

validEstMethods = {'QR','fitlm'};
if ~ismember(cfg.estimationMethod, validEstMethods)
    error('vba_cfg: cfg.estimationMethod ''%s'' is not valid. Choose from: %s.', ...
        cfg.estimationMethod, strjoin(validEstMethods, ', '));
end

if cfg.doRobust && strcmp(cfg.estimationMethod, 'QR')
    warning('vba_cfg: doRobust=1 requires estimationMethod=''fitlm''. Switching automatically.');
    cfg.estimationMethod = 'fitlm';
end

validFitMethods = {'REML','ML'};
if ~ismember(cfg.lme_fitMethod, validFitMethods)
    error('vba_cfg: cfg.lme_fitMethod ''%s'' is not valid. Use ''REML'' or ''ML''.', cfg.fitMethod);
end

% =========================================================================
% Print summary if running interactively
% =========================================================================
if ~cfg.doSlurm
    vba_cfg_print(cfg);
end

end


% -------------------------------------------------------------------------
% Local helpers
% -------------------------------------------------------------------------

function cfg = default(cfg, field, val)
% Set cfg.field = val only if field is not already present
if ~isfield(cfg, field) || isempty(cfg.(field))
    cfg.(field) = val;
end
end

function cfg = require(cfg, field, msg)
% Error if a required field is absent
if ~isfield(cfg, field) || isempty(cfg.(field))
    error('vba_cfg: %s', msg);
end
end

function vba_cfg_print(cfg)
% Print a readable summary of the active configuration
fprintf('\n--- VBA Configuration ---\n');
fprintf('  Model      : %s\n', cfg.model);
fprintf('  Model type : %s\n', cfg.modelType);
fprintf('  Estimation : %s\n', cfg.estimationMethod);
fprintf('  Mask       : %s\n', cfg.f_mask);
fprintf('  Permutations: %d (starting at %d)\n', cfg.numPerm, cfg.startPerm);
fprintf('  Z-score    : %d  |  Log-transform: %d  |  Robust: %d\n', ...
    cfg.doZscore, cfg.doLogTrans, cfg.doRobust);
if strcmp(cfg.modelType, 'mixed')
    fprintf('  Random effects: %s  |  Fit method: %s\n', ...
        cfg.randomEffects, cfg.fitMethod);
end
fprintf('-------------------------\n\n');
end