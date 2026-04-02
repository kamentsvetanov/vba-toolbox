function [varout] = vba_prepare_model(Model)
% Parse Wilkinson formula, identify covariates, numeric/categorical variables
% This code extracts variable names from an equation defined using
% Wilkinsons' notation


% -------------------------
% Get Variable names
% -------------------------
VarNames = strtrim(split(Model,["~","*","+","^",":"]));
ResponseVar = VarNames(1); 
PredictorVars = VarNames(2:end);


% -------------------------
% Create output file where Response and Predictors are prefixed by 'dv' and
% 'iv', respectively.
% -------------------------
strModel = regexprep(Model,'f_','');
strModel = regexprep(strModel,'c_','');
strModel = regexprep(strModel,'(\<[a-z])','${upper($1)}');
strModel = regexprep(strModel,{'+',' '},'');
strModel = regexprep(strModel,'~','_iv');
strModel = regexprep(strModel,'*','X');
strModel = ['dv' strModel];
% strModel = ['dv' strModel '_' datestr(now,'yyyymmdd')];


% ----------------------------
% Check for interaction terms
% (Works with 1 interaction term at the moment)
% ----------------------------
if contains(Model,{'*',':'})
    moderation.Y = ResponseVar{1};
    
    % Check for 3 way interaction
    str = regexp(Model,'(\w*)[*|:](\w*)[*|:](\w*)','tokens');
    
    % Check for 2-way interaction
    if ~isempty(str)
        moderation.X = str{1}{1};
        moderation.W = str{1}{2};
        moderation.Z = str{1}{3};
        moderation.C = PredictorVars(~ismember(PredictorVars,[moderation.X,moderation.W,moderation.Z]));
    else
        str = regexp(Model,'(\w*)[*|:](\w*)','tokens');
        
        if ~isempty(str)
            moderation.X = str{1}{1};
            moderation.Z = str{1}{2};
            moderation.C = PredictorVars(~ismember(PredictorVars,[moderation.X,moderation.Z]));
        end
    end
        
%     moderation.X = char(string(regexp(Model,'(\w*)[*?:?](?!**)','tokens'))); % Predictor is before '*' or ':' symbol
%     moderation.Z = char(string(regexp(Model,'[*?:?](\w*)','tokens'))); % Moderator is before '*' or ':' symbol
    
    varout.moderation = moderation;
end

% -------------------------
% Prepare Output structure
varout.ResponseVar      = ResponseVar;
varout.PredictorVars    = PredictorVars;
varout.f_model          = strModel; 
varout.VarNames         = VarNames;