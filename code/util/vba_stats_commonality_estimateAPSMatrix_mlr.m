function APSMatrix = vba_stats_commonality_estimateAPSMatrix_mlr(dataTemp,Model,numcc,doRobust)

APSMatrix = nan(numcc,1);
for i = 1:numcc
%         tic
    mlrTemp  = fitlm(dataTemp,Model{i},'RobustOpts',doRobust);
%         toc
    APSMatrix(i)    = mlrTemp.Rsquared.Ordinary;
%         APSMatrix{i, 'R2'} = mlrTemp.Rsquared.Ordinary;
%         APSMatrix{i, 'k'} = sum(PredBitMap{:,i});
    
    %-(not implemented) Possibly a more efficient way to estimate Total variance explained by the
    % model. Note that it does not work for interactions or squared
    % terms defined in Wilkinson annotation. Instead these should be
    % modelled by the user.
    
%         y = dataTemp.(dv);
%         X = [ones(Ns,1) dataTemp{:,ModelPred{i}}];
%         mu_y    = mean(y);
%         if doRobust
%             [betas stats]= robustfit(X,y,[],[],'off');
%             ss_res  = sum(stats.resid.^2);
%         else
% %             [b,bint,r,rint,stats] = regress(dataTemp.(dv),dataTemp{:,ModelPred{i}});
%             betas   = X\y; % add constant term to make identical to robust fit and fitlm
%             ss_res  = sum((X*betas - y).^2);
%         end
%         ss_tot  = sum((y - mu_y).^2);
%         R2      = 1 - ss_res/ss_tot;
%         APSMatrix(i) = 1 - ss_res/ss_tot;
%         tic
%         y = dataTemp.(dv);
%         X = dataTemp{:,ModelPred{i}};
%         mu_y    = mean(y);
%         betas   = x\y;
%         ss_tot  = sum((y - mu_y).^2);
%         ss_res  = sum((X*betas - y).^2);
%         R2      = 1 - ss_res/ss_tot;
%         toc
    
    
end