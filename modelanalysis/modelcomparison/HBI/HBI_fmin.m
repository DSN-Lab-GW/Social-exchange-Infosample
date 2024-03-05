fname = load('lap_Threshold.mat');
cbm   = fname.cbm;
% look at fitted parameters
cbm.output.parameters

fname = load('lap_Optimal.mat');
cbm   = fname.cbm;
% look at fitted parameters
cbm.output.parameters

fname = load('lap_Uncertainty.mat');
cbm   = fname.cbm;
% look at fitted parameters
cbm.output.parameters
%% 
load lap_Threshold.mat
load estimates_DDM_TD_all pars_est
for i = 1:size(cbm.output.parameters,1)
    to_mean(i,1) = 1/(1+exp(-cbm.output.parameters(i,1)));
    to_mean(i,2) = 25/(1+exp(-cbm.output.parameters(i,2)));
    to_mean(i,3) = 25/(1+exp(-cbm.output.parameters(i,3)));
end
for i = 1:3
   [rho(i),pval(i)] = corr(pars_est(:,i+1),to_mean(:,i));
   diff(:,i) = pars_est(:,i+1) - to_mean(:,i);
   [h,p(i),ci,stats] = ttest(diff(:,i));
end
% threshold
fig = figure(6)
i=1
subplot(1,size(to_mean,2),i)
xlabel('Input parameter')
ylabel('Recovered parameter')
scatter(pars_est(:,i+1),to_mean(:,i))
% xlim([0 4]);
% ylim([0 4]);
xlabel('Decision noise')
hold on
i=2
subplot(1,size(to_mean,2),i)
scatter(pars_est(:,i+1),to_mean(:,i))
% xlim([0 100]);
% ylim([0 100]);
xlabel('Positive bound')
i=3
subplot(1,size(to_mean,2),i)
scatter(pars_est(:,i+1),to_mean(:,i))
% xlim([0 100]);
% ylim([0 100]);
xlabel('Negative bound')
% subplot label/title
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'HBI');
title(han,'HBI vs fmin');
hold off


%Uncertainty
% load Uncertainty_recov_1.mat Uncertainty_Uncertainty
% to_mean = Uncertainty_Uncertainty.pars_est;
% load Uncertainty_recov_2.mat Uncertainty_Uncertainty
% to_mean = (to_mean + Uncertainty_Uncertainty.pars_est)/2;
load lap_Uncertainty.mat
load estimates_uncertainty_TD_all pars_est
for i = 1:size(cbm.output.parameters,1)
    to_mean(i,1) = exp(cbm.output.parameters(i,1));
    to_mean(i,2) = 1/(1+exp(-cbm.output.parameters(i,2)));
    to_mean(i,3) = 1+exp(cbm.output.parameters(i,3));
    to_mean(i,4) = 1+exp(cbm.output.parameters(i,4));
end
for i = 1:4
   [rho(i),pval(i)] = corr(pars_est(:,i+1),to_mean(:,i));
   diff(:,i) = pars_est(:,i+1) - to_mean(:,i);
   [h,p(i),ci,stats] = ttest(diff(:,i));
end

fig = figure(7)
i=1
subplot(1,size(to_mean,2),i)
xlabel('Input parameter')
ylabel('Recovered parameter')
scatter(pars_est(:,i+1),to_mean(:,i))
xlim([0 500]);
ylim([0 500]);
xlabel('Decision noise')
hold on
i=2
subplot(1,size(to_mean,2),i)
scatter(pars_est(:,i+1),to_mean(:,i))
xlim([0 0.15]);
ylim([0 0.15]);
xlabel('Uncertainty criterion')
i=3
subplot(1,size(to_mean,2),i)
scatter(pars_est(:,i+1),to_mean(:,i))
xlim([0 20]);
ylim([0 20]);
xlabel('Alphaprior')
i=4
subplot(1,size(to_mean,2),i)
scatter(pars_est(:,i+1),to_mean(:,i))
xlim([0 20]);
ylim([0 20]);
xlabel('Betaprior')
% subplot label/title
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'HBI');
title(han,'HBI vs fmin');
hold off

%Optimal
% load Optimal_recov_1.mat Optimal_Optimal
% to_mean = Optimal_Optimal.pars_est;
% load Optimal_recov_2.mat Optimal_Optimal
% to_mean = (to_mean + Optimal_Optimal.pars_est)/2;
load lap_Optimal.mat
load estimates_optimal_TD_all pars_est
for i = 1:size(cbm.output.parameters,1)
    to_mean(i,1) = -1/(1+exp(-cbm.output.parameters(i,1)));
    to_mean(i,2) = 10+100/(1+exp(-cbm.output.parameters(i,2)));
    to_mean(i,3) = 5/(1+exp(-cbm.output.parameters(i,3)));
    to_mean(i,4) = (1/(1+exp(-cbm.output.parameters(i,4)))-0.5)*2;
    to_mean(i,5) = 1+25/(1+exp(-cbm.output.parameters(i,5)));
    to_mean(i,6) = 1+25/(1+exp(-cbm.output.parameters(i,6)));
end
for i = 1:6
   [rho(i),pval(i)] = corr(pars_est(:,i+1),to_mean(:,i));
   diff(:,i) = pars_est(:,i+1) - to_mean(:,i);
   [h,p(i),ci,stats] = ttest(diff(:,i));
   subplot(1,size(to_mean,2),i)
    plot(diff(:,i))
end

fig = figure(8)
i=1
subplot(1,size(to_mean,2),i)
ylabel('Recovered parameter')
scatter(pars_est(:,i+1),to_mean(:,i))
xlim([-0.01 0.01]);
ylim([-0.01 0.01]);
xlabel('k')

i=2
subplot(1,size(to_mean,2),i)
scatter(pars_est(:,i+1),to_mean(:,i))
xlim([0 10000]);
ylim([0 10000]);
xlabel('Decision noise')
i=3
subplot(1,size(to_mean,2),i)
scatter(pars_est(:,i+1),to_mean(:,i))
xlim([0 1]);
ylim([0 1]);
xlabel('Cost')
hold on
i=4
subplot(1,size(to_mean,2),i)
scatter(pars_est(:,i+1),to_mean(:,i))
xlim([0 1]);
ylim([0 1]);
xlabel('Risk')
i=5
subplot(1,size(to_mean,2),i)
scatter(pars_est(:,i+1),to_mean(:,i))
xlim([0 25]);
ylim([0 25]);
xlabel('Alphaprior')
i=6
subplot(1,size(to_mean,2),i)
scatter(pars_est(:,i+1),to_mean(:,i))
xlim([0 25]);
ylim([0 25]);
xlabel('Betaprior')
% subplot label/title
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'HBI');
title(han,'HBI vs fmin');
hold off

%%
load lap_Optimal.mat
load estimates_optimal_TD_all pars_est
for i = 1:size(cbm.output.parameters,1)
    to_mean(i,1) = -1/(1+exp(-cbm.output.parameters(i,1)));
    to_mean(i,2) = 10+100/(1+exp(-cbm.output.parameters(i,2)));
    to_mean(i,3) = 5/(1+exp(-cbm.output.parameters(i,3)));
    to_mean(i,4) = (1/(1+exp(-cbm.output.parameters(i,4)))-0.5)*2;
    to_mean(i,5) = 1+25/(1+exp(-cbm.output.parameters(i,5)));
    to_mean(i,6) = 1+25/(1+exp(-cbm.output.parameters(i,6)));
end

badfit_num = [1 5 7 8 9 16 25 26 28 30];
allfit_num = [1:36];
goodfit_num = allfit_num(setdiff(1:end,badfit_num));
for i = 1:length(badfit_num)
    badfit(i,:) = to_mean(badfit_num(i),:);
end
for i = 1:length(goodfit_num)
    goodfit(i,:) = to_mean(goodfit_num(i),:);
end

fig = figure(1)
subplot(2,1,1)
boxplot(goodfit)
ylabel('Good fit');

subplot(2,1,2)
boxplot(badfit)
ylabel('Bad fit');
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
xlabel(han,'Optimal model parameter');
title(han,'Good fit vs Bad fit');

%%
load lap_Optimal.mat
load estimates_optimal_TD_all pars_est
for i = 1:size(cbm.output.parameters,1)
    to_mean(i,1) = -1/(1+exp(-cbm.output.parameters(i,1)));
    to_mean(i,2) = 10+100/(1+exp(-cbm.output.parameters(i,2)));
    to_mean(i,3) = 5/(1+exp(-cbm.output.parameters(i,3)));
    to_mean(i,4) = (1/(1+exp(-cbm.output.parameters(i,4)))-0.5)*2;
    to_mean(i,5) = 1+25/(1+exp(-cbm.output.parameters(i,5)));
    to_mean(i,6) = 1+25/(1+exp(-cbm.output.parameters(i,6)));
end

fig = figure(1)
subplot(2,1,1)
boxplot(to_mean)
ylabel('Good fit');

subplot(2,1,2)
boxplot(pars_est(2:end))
ylabel('Bad fit');
han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
xlabel(han,'Optimal model parameter');
title(han,'Good fit vs Bad fit');
