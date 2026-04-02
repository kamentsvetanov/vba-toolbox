function vba_write_results(cfg, tvals, nameOutput, Vdv, idxMask, iperm)

% tvals, Vdv, idxMask, nameOutput, outDir, iperm)
outDir = cfg.outDir;

for i=1:numel(nameOutput)
%                 % Output coefficients not containing covariates 
%                 % (i.e. predictors prefixed by 'c_')
   % Do not create folders for for a pattern specified in excludeFromResult
    if ~contains(nameOutput{i},cfg.excludeFromResults)
        
        
        Vtemp           = Vdv;
        Vtemp.pinfo(1)  = 1;
        tmap            = zeros(Vtemp.dim);
%         bmap            = zeros(Vtemp.dim);

        % Save t-stats
        outname         = regexprep(nameOutput{i},'f_','');
        tdir            = fullfile(outDir,['tval_',outname]); 
        tmap(idxMask)   = tvals(:,i);
       
        Vtemp.fname     = fullfile(tdir,sprintf('results_null_%.5d.nii',iperm));
        
        spm_write_vol(Vtemp,tmap);

        % Save beta coefficients (r-values if data was z-scoredd)
%         bdir            = fullfile(outDir,['bval_',nameOutput]); 
%         bmap(idxMask)   = bvals(:,i);  
%         Vtemp.fname     = fullfile(bdir,sprintf('results_null_%.5d.nii',iperm));
%         spm_write_vol(Vtemp,bmap);
    end
  
end