% FITTING for ORA study
clear all; close all
% %data
% load data_ThresholdC_sim_ili.mat
% cd /lustre/groups/rosenblaugrp/4HPC_infosample/model/modelrecovery_code
load data_ThresholdC_sim.mat
%% fit uncertainty model
model_num = 111
subjvec = size(data_generate,2);
numinit = 100; % Number of starting points for parameter fitting
for num_ite = 1:3%size(data_generate,1)
    num_ite
    for subjidx = 1:subjvec 
%         subjidx
        datasubj.red    = data_generate{num_ite,subjidx}.red;
        datasubj.green  = data_generate{num_ite,subjidx}.green;
        datasubj.choice = data_generate{num_ite,subjidx}.choice;

        myNLL = @(pars) mymodelUncertainty_ORA(pars, datasubj);

        init         = NaN(numinit, 2); % nr of parameters 
        init(:,1)    = log(randi(150, numinit,1)); %beta softmax
        init(:,2)    = rand(numinit,1); % criterion
        init(:,3)     = rand(numinit,1) + 1; %Alphaprior
        init(:,4)     = rand(numinit,1) + 1; %Betaprior

        lowLimits =  [0   0 1  1 ];
        highLimits = [inf 1 inf inf];

        for runidx = 1:numinit
            [pars_per_run(subjidx, runidx, :), NLL(runidx)] = fmincon(myNLL, init(runidx,:),[],[],[],[], lowLimits, highLimits, [], optimset('Display', 'off'));
        end
%         NLL
        [~, bestrun] = min(NLL); 
        [fittedpars, bestNLL] = fmincon(myNLL, init(bestrun,:),[],[],[],[], lowLimits, highLimits, [], optimset('Display', 'off'));

        pars_est_uncertainty{num_ite}(subjidx,:) = fittedpars;
        allbestNLL_uncertainty{num_ite}(subjidx) = bestNLL';

    end
end
% save('/CCAS/home/wenda_liu/threC_uncer.mat','pars_est_uncertainty','allbestNLL_uncertainty')
save('threC_uncer.mat','pars_est_uncertainty','allbestNLL_uncertainty');
clearvars -except data_generate subjvec numinit pars_est_uncertainty allbestNLL_uncertainty
%% fit thresholdC model
model_num = 222
subjvec = size(data_generate,2);
numinit = 100; % Number of starting points for parameter fitting
for num_ite = 1:3%size(data_generate,1)
    num_ite
    for subjidx = 1:subjvec 
        subjidx
        datasubj.red    = data_generate{num_ite,subjidx}.red;
        datasubj.green  = data_generate{num_ite,subjidx}.green;
        datasubj.choice = data_generate{num_ite,subjidx}.choice;

        myNLL = @(pars) mymodelDDM_twobounds_ORA(pars, datasubj);

        init         = NaN(numinit, 3);
        init(:,1)    = rand(numinit,1); %beta softmax
        init(:,2)    = randi(13, numinit,1); % bound
        init(:,3)    = -rand(numinit,1); % collapse rate

        lowLimits = [0 0 -1];
        highLimits = [100 inf 0];

         for runidx = 1:numinit
        [pars_per_run(subjidx, runidx, :), NLL(runidx)] = fmincon(myNLL, init(runidx,:),[],[],[],[], lowLimits, highLimits, [], optimset('Display', 'on'));
        end
%         NLL
        [~, bestrun] = min(NLL); 
        [fittedpars, bestNLL] = fmincon(myNLL, init(bestrun,:),[],[],[],[], lowLimits, highLimits, [], optimset('Display', 'off'));

        pars_est_thresholdC{num_ite}(subjidx,:) = fittedpars;
        allbestNLL_thresholdC{num_ite}(subjidx) = bestNLL';

    end
end
% save('/CCAS/home/wenda_liu/threC_threC.mat','pars_est_thresholdC','allbestNLL_thresholdC')
save('threC_threC.mat','pars_est_thresholdC','allbestNLL_thresholdC');
clearvars -except data_generate subjvec numinit pars_est_uncertainty allbestNLL_uncertainty pars_est_thresholdC allbestNLL_thresholdC
%% fit optimal model
model_num = 333
subjvec = size(data_generate,2);
numinit = 100; % Number of starting points for parameter fitting
for num_ite = 1:3%size(data_generate,1)
    num_ite
    for subjidx = 1:subjvec
%         subjidx
        datasubj.red    = data_generate{num_ite,subjidx}.red;
        datasubj.green  = data_generate{num_ite,subjidx}.green;
        datasubj.choice = data_generate{num_ite,subjidx}.choice;

        myNLL = @(pars) mymodel_Optimal_ORA(pars, datasubj);

        init        = NaN(numinit, 6); 
        init(:,1)   = rand(numinit,1) * -0.1;      %k
        init(:,2)   = randi(150, numinit,1) + 50;  %beta 0
        init(:,3)   = rand(numinit,1) * 0.1;       %cost
        init(:,4)   = rand(numinit,1) * 0.1;       %risk
        init(:,5)   = randn(numinit,1) + 1;        %Alphaprior
        init(:,6)   = randn(numinit,1) + 1;        %Betaprior

    %     lowLimits =  [-1  10  0];
    %     highLimits = [ 0  inf 10]; 
    %     
    %     lowLimits =  [-1  10  0  -1];
    %     highLimits = [ 0  inf 10   1]; 
    %     
        lowLimits =  [-1  10  0  -1  1   1];
        highLimits = [ 0  inf 5   1  26 26]; 

        for runidx = 1:numinit
            [pars_per_run(subjidx, runidx, :), NLL(runidx)] = fmincon(myNLL, init(runidx,:),[],[],[],[], lowLimits, highLimits, [], optimset('Display', 'off'));
        end
%         NLL
        [~, bestrun] = min(NLL);
        [fittedpars, bestNLL] = fmincon(myNLL, init(bestrun,:),[],[],[],[], lowLimits, highLimits, [], optimset('Display', 'off'));


        pars_est_optimal{num_ite}(subjidx,:) = fittedpars;
        allbestNLL_optimal{num_ite}(subjidx) = bestNLL';

    end
end
% save('/CCAS/home/wenda_liu/threC_optimal.mat','pars_est_optimal','allbestNLL_optimal')
save('threC_optimal.mat','pars_est_optimal','allbestNLL_optimal');
clearvars -except data_generate subjvec numinit pars_est_uncertainty allbestNLL_uncertainty pars_est_thresholdC allbestNLL_thresholdC pars_est_optimal allbestNLL_optimal
%%
for num_ite = 1:size(data_generate,1)
    for subjidx = 1:subjvec
        dpssym{num_ite}(subjidx,1) = subjidx;
        dpssym{num_ite}(subjidx,2) = size(data_generate{num_ite,subjidx}.green,1);
    end
    %load model estimates
    uncertNLL = allbestNLL_uncertainty{num_ite}'; 
    thresholdNLL =allbestNLL_thresholdC{num_ite}';
    optimalNLL = allbestNLL_optimal{num_ite}';
    %AIC's & BIC's
    optimalAIC  = -2* (( -optimalNLL) - 6);
    optimalBIC  = (log(dpssym{num_ite}(:,2)))*6 - (2* -optimalNLL);

    uncertAIC  = -2* (( -uncertNLL) - 4); 
    uncertBIC  = (log(dpssym{num_ite}(:,2)))*4 - (2* -uncertNLL);

    countAIC  = -2* (( -thresholdNLL) - 3);
    countBIC  = (log(dpssym{num_ite}(:,2)))*3 - (2* -thresholdNLL);
    %compare sum BIC for each iteration
    All_AIC = [uncertAIC countAIC optimalAIC];
    All_BIC = [uncertBIC countBIC optimalBIC];
    [L_value,L_index(num_ite)] = min(sum(All_BIC,1)'); %uncertainty 1;optimal  3; threshold 2
end

freq_win(1) = numel(find(L_index == 1));%uncertainty 1;
freq_win(2) = numel(find(L_index == 2));%threshold 2
freq_win(3) = numel(find(L_index == 3));%optimal 3;

bar(freq_win)
names = {'Uncertainty','Threshold','Optimal'};
set(gca,'xticklabel',names)
ylabel('Win frequency')
% save('/CCAS/home/wenda_liu/threC_sim_recov_BIC.mat','freq_win','All_BIC','All_AIC')
save('threC_sim_recov_BIC.mat','freq_win','All_BIC','All_AIC')
%%
% clear all
% load optimal_optimal.mat
% for i = 1:size(pars_est_optimal{1},2)
%     for j = 1:size(pars_est_optimal,2)
%         abc{i}(:,j) = pars_est_optimal{j}(:,i);
%     end
% end
% clear all
% load threC_threC.mat
% for i = 1:size(pars_est_thresholdC{1},2)
%     for j = 1:size(pars_est_thresholdC,2)
%         abc{i}(:,j) = pars_est_thresholdC{j}(:,i);
%     end
% end
% clear all
% load uncer_uncer.mat
% for i = 1:size(pars_est_uncertainty{1},2)
%     for j = 1:size(pars_est_uncertainty,2)
%         abc{i}(:,j) = pars_est_uncertainty{j}(:,i);
%     end
% end
% %
% for i = 1:size(abc,2)
%     to_plot{i}(:,1) = mean(abc{i}(:,1:50),2);
%     to_plot{i}(:,2) = mean(abc{i}(:,51:100),2);
% end
% for i = 1:size(abc,2)
%     subplot(1,size(to_plot,2),i)
%     scatter(to_plot{i}(:,1),to_plot{i}(:,2))
% end
% for i = 1:size(abc,2)
% % [b,bint,r,rint,stats] = regress(to_plot{i}(:,1),[ones(size(to_plot{i}(:,2))) to_plot{i}(:,2)]);
% mdl = fitlm(to_plot{i}(:,1),to_plot{i}(:,2))
% end
% mdl = fitlm(to_plot{i}(:,1),to_plot{i}(:,2))