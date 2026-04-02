function [Y, meta] = vba_load_images(tbl, VarNames, cfg)
%
% VBA_LOAD_IMAGES  Load and prepare imaging data for voxel-based analysis.
%
% Loads brain mask, imports all imaging modalities specified in VarNames
% (identified by 'f_' prefix), applies optional log-transform, and
% computes a valid-voxel mask across all subjects and modalities.
%
% INPUTS
%   tbl       - Subject table containing image filepaths and covariates
%   VarNames  - Cell array of variable names from the model (parsed)
%   cfg.f_mask- Filepath to brain mask NIfTI image
%   cfg       - Configuration struct with optional fields:
%                 cfg.doLogTrans  (default 0)
%
% OUTPUTS
%   Y         - [numSub x numValidVox x numMap] imaging data array,
%               restricted to voxels valid across all subjects and maps
%   meta      - Struct containing:
%                 meta.Ymask        - [x y z] logical brain mask volume
%                 meta.dimMask      - dimensions of mask [x y z]
%                 meta.idxMask      - linear indices of mask voxels in 3D volume
%                 meta.idxValidVox  - indices into idxMask of valid voxels
%                 meta.idxMaskValid - linear indices into 3D volume (composed)
%                 meta.numVox       - number of valid voxels (post NaN filter)
%                 meta.numSub       - number of subjects
%                 meta.numMap       - number of imaging modalities
%                 meta.VarNamesMaps - cell array of imaging variable names
%                 meta.VarNamesGlobal - cell array of non-imaging variable names
%                 meta.Vdv          - SPM volume header of first/DV modality
%
% See also: vba_run, spm_vol, spm_read_vols
%
% Author : Kamen Tsvetanov, Ph.D.
% Affil. : Department of Clinical Neurosciences, University of Cambridge
%__________________________________________________________________________

% -------------------------
% Parse config
% -------------------------
try doLogTrans = cfg.doLogTrans; catch, doLogTrans = 0; end % Whether or not to log-transform neuroimaging data

% -----------------------------------------------
% Identify imaging vs global variable names
% -----------------------------------------------
idxMaps         = find(contains(VarNames, 'f_'));
VarNamesMaps    = VarNames(idxMaps);
VarNamesGlobal  = VarNames(~contains(VarNames, 'f_'));

numMap  = numel(idxMaps);
numSub  = size(tbl, 1);

if numMap == 0
    error('vba_load_images: No imaging variables (f_ prefix) found in VarNames.');
end

% -----------------------------------------------
% Load brain mask
% -----------------------------------------------
fprintf('Loading brain mask: %s\n', f_mask);
Vmask   = spm_vol(f_mask);
Ymask   = logical(spm_read_vols(Vmask));
dimMask = size(Ymask);
idxMask = find(Ymask);
numVoxMask = numel(idxMask);

fprintf('  Mask contains %d voxels\n', numVoxMask);

% -----------------------------------------------
% Load each imaging modality
% -----------------------------------------------
% Pre-allocate over full mask first; valid-voxel
% filtering is applied after all maps are loaded
Y_full = nan(numSub, numVoxMask, numMap);
Vdv    = [];

for iMap = 1:numMap

    nameMaps = VarNames{idxMaps(iMap)};
    fprintf('Loading modality %d/%d: %s\n', iMap, numMap, nameMaps);

    % Resolve filepaths from table
    fpaths = tbl.(nameMaps);
    if iscell(fpaths)
        fpaths = char(fpaths);
    end

    V = spm_vol(fpaths);

    % Dimension check against mask
    if ~isequal(dimMask, V(1).dim)
        error('vba_load_images: Dimension mismatch between mask [%s] and %s [%s].', ...
            num2str(dimMask), nameMaps, num2str(V(1).dim));
    end

    % Read volumes and reshape to [numSub x numVoxMask]
    y = spm_read_vols(V);           % [x y z numSub]
    y = permute(y, [4 1 2 3]);      % [numSub x y z]
    Y_full(:, :, iMap) = y(:, Ymask);

    % Store header of first map (assumed to be DV or reference space)
    if iMap == 1
        Vdv = V(1);
    end
end

% -----------------------------------------------
% Optional log-transform (applied before NaN check)
% -----------------------------------------------
% N.B. Applied to all modalities. Negative or zero values will produce
% NaN/Inf ? these will be caught by the valid-voxel filter below.
if doLogTrans
    fprintf('Applying log-transform to imaging data\n');
    Y_full = log(Y_full);
end

% -----------------------------------------------
% V