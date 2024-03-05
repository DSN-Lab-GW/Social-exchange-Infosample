
%model comparisons
load('summaryORA3.mat') %load datafile
dpssym = decisionsPerSubject(data, []);


%before running, check if 2nd column is indeed the NLL in your data
Uncert = load('estimates_uncertaintyFull3_ORA'); %load Uncertainty model estimates
uncertNLL = Uncert.allbestNLL(:,1); 

DDM = load('estimates_DDMcollapse_ORA3.mat'); %load threshold model estimates
DDMNLL = DDM.allbestNLL(:,1);

optimal = load('estimates_optimalFull_ORA3'); %load sample cost model estimates
optimalNLL = optimal.allbestNLL(:,1);

total = load('estimates_TOTAL.mat');
totalNLL = total.allbestNLL(:,1);
%
% DDM2C = load('estimates_DDMTwoboundsCollapse_ORA3_20210308.mat'); %load threshold model estimates
% DDM2CNLL = DDM2C.allbestNLL(:,1);
% 
% DDM2 = load('estimates_DDMTwobounds_ORA3_20210308.mat'); %load sample cost model estimates
% DDM2NLL = DDM2.allbestNLL(:,1);
% 
% DDMC = load('estimates_DDMcollapse_ORA3_20210308'); %load Uncertainty model estimates
% DDMCNLL = DDMC.allbestNLL(:,1); 
% 
% DDMB = load('estimates_DDMbasic_ORA3_20210308'); %load Uncertainty model estimates
% DDMBNLL = DDMB.allbestNLL(:,1); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % AIC's & BIC's

 % THIS IS BETWEEN MODELS. First do the within model comparisons to see
 % which is the best model. The code will need to be adjusted but to load
 % the correct files
 
optimalAIC  = -2* (( -optimalNLL) - 6);
optimalBIC  = (log(dpssym(:,2)))*6 - (2* -optimalNLL);

uncertAIC  = -2* (( -uncertNLL) - 4); 
uncertBIC  = (log(dpssym(:,2)))*4 - (2* -uncertNLL);

countAIC  = -2* (( -DDMNLL) - 3);
countBIC  = (log(dpssym(:,2)))*3 - (2* -DDMNLL);

totalAIC  = -2* (( -totalNLL) - 2);
totalBIC  = (log(dpssym(:,2)))*2 - (2* -totalNLL);
%
% countBAIC  = -2* (( -DDMBNLL) - 2);
% countBBIC  = (log(dpssym(:,2)))*2 - (2* -DDMBNLL);
% 
% countCAIC  = -2* (( -DDMCNLL) - 3);
% countCBIC  = (log(dpssym(:,2)))*3 - (2* -DDMCNLL);
% 
% count2AIC  = -2* (( -DDM2NLL) - 3);
% count2BIC  = (log(dpssym(:,2)))*3 - (2* -DDM2NLL);
% 
% count2CAIC  = -2* (( -DDM2CNLL) - 4);
% count2CBIC  = (log(dpssym(:,2)))*4 - (2* -DDM2CNLL);

% bayesianValue = ( numTrials * log(SSE/numTrials) ) + ( numParam * log(numTrials) );
%here's where you calculate the difference in AIC or BICs between models
diffsA =  uncertAIC - optimalAIC; 
diffsB =  uncertBIC - optimalBIC; 
% diffsA =  optimalAIC - countAIC; 
% diffsB =  optimalBIC - countBIC; 

% AIC bootstrapped 95% CI of the summed difference
AIC = diffsA; 
nboot = 10000000;
sum_boot = NaN(nboot,1);
for i = 1:nboot
   % Draw bootstrap sample
   AIC_boot = randsample(AIC, length(AIC), 1);
   sum_boot(i) = sum(AIC_boot);
end
sumAIC = sum(AIC)
CI_lowerAIC = quantile(sum_boot, 0.025);
CI_upperAIC = quantile(sum_boot, 0.975);
disp('CI_Low_AIC')
disp(num2str(CI_lowerAIC,'%.2f'))
disp('CI_Up_AIC')
disp(num2str(CI_upperAIC,'%.2f'))

% BIC bootstrapped 95% CI of the summed difference
BIC = diffsB; 
nboot = 1000000;
sum_boot = NaN(nboot,1);
for i = 1:nboot
   % Draw bootstrap sample
   BIC_boot = randsample(BIC, length(BIC), 1);
   sum_boot(i) = sum(BIC_boot);
end
sumBIC = sum(BIC)
CI_lowerBIC = quantile(sum_boot, 0.025);
CI_upperBIC = quantile(sum_boot, 0.975);
disp('CI_Low_BIC')
disp(num2str(CI_lowerBIC,'%.2f'))
disp('CI_Up_BIC')
disp(num2str(CI_upperBIC,'%.2f'))

%%
All_AIC = [optimalAIC uncertAIC countAIC];
All_BIC = [optimalBIC-countBIC uncertBIC-countBIC countBIC-countBIC totalBIC-countBIC];
figure(1)
% bar(All_AIC,'DisplayName','to_plot')%/All_BIC
b = barh(sum(All_BIC,1),'FaceColor','flat')
b.CData(1,:) = [0.5 0.5 0.5];
b.CData(2,:) = [0.7 0.7 0.7];
b.CData(3,:) = [0.9 0.9 0.9];
% xlabel('Sample cost','Uncertainty','Threshold')
xlabel('\Delta BIC')%/Sum BIC
% legend('Sample cost','Uncertainty','Threshold')
% ylim([0 500])
yticks([1 2 3])
yticklabels({'Sample cost','Uncertainty','Threshold'})

%% model position
num_model_position(1,1) = length(find((optimalBIC < uncertBIC)&(optimalBIC < countBIC)))
num_model_position(2,1) = length(find((uncertBIC < optimalBIC)&(uncertBIC < countBIC)))
num_model_position(3,1) = length(find((countBIC < uncertBIC)&(countBIC < optimalBIC)))

num_model_position(1,3) = length(find((optimalBIC > uncertBIC)&(optimalBIC > countBIC)))
num_model_position(2,3) = length(find((uncertBIC > optimalBIC)&(uncertBIC > countBIC)))
num_model_position(3,3) = length(find((countBIC > uncertBIC)&(countBIC > optimalBIC)))

num_model_position(:,2) = numel(sub_code);
num_model_position(:,2) = num_model_position(:,2) - num_model_position(:,1) - num_model_position(:,3)

figure(2)
bb = bar(num_model_position')
bb.CData(1,1) = [0.5 0.5 0.5];
bb.CData(2,:) = [0.3 0.3 0.3];
bb.CData(3,:) = [0.1 0.1 0.1];
ylim([0 numel(sub_code)]);
names = {'First','Second','Third'};
set(gca,'xticklabel',names)
legend('Sample cost','Uncertainty','Threshold','Location','Northeast')
ylabel('Num of subject')
xlabel('')

% %% BIC function
% function bayesianValue = calcBIC(numTrials, numParam, SSE)
%     bayesianValue = ( numTrials * log(SSE/numTrials) ) + ( numParam * log(numTrials) );
% end
%% BMS
i_sit = 1;
[BMS.alpha{ i_sit }, BMS.exp_r{ i_sit }, BMS.xp{ i_sit }, BMS.pxp{ i_sit }, BMS.bor{ i_sit }] = spm_BMS( -All_BIC(:,[1,2,3]) ); 
figure(3)
hBar = bar(BMS.pxp{1,1});
set(hBar, 'FaceColor','0 0.4470 0.7410')
hold on
names = {'Sample cost','Uncertainty','Threshold'};
set(gca,'xticklabel',names)
ylabel('Exceedance probability')
ylim([0 1]);
hold off