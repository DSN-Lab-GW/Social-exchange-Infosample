%% prepare data
close all;clear all
load summaryORA3_ASD_male_matched.mat
load summaryORA3_TD.mat

data_all = data;
clear data
subjvec = unique(data_all.subjid);
for subjidx = 1:length(subjvec) % for generated data run 20 subjects
    subjidx
    idx = find(data_all.subjid == subjvec(subjidx));
    data{subjidx,1}.red    = data_all.red(idx);
    data{subjidx,1}.green  = data_all.green(idx);
    data{subjidx,1}.choice = data_all.choice(idx);
end
%% check individual data
% subj1 = data{3};
% parameters = randn(1,3);
% parameters = [1 10 10]
% F1 = mymodelDDM_twobounds_ORA(parameters,subj1)
% aa = repmat(100,1,100)
% bb = aa(1)/(1+exp(-parameters(1)))
%% set parameter initiation
v     = 6.25;
% prior_Threshold = struct('mean',zeros(3,1),'variance',v); % note dimension of 'mean'
prior_Optimal = struct('mean',zeros(4,1),'variance',v); % note dimension of 'mean'
prior_Uncertainty = struct('mean',zeros(4,1),'variance',v); % note dimension of 'mean'
%% output file
% fname_Threshold = 'lap_Threshold.mat'; 
fname_Optimal = 'lap_Optimal.mat'; 
fname_Uncertainty = 'lap_Uncertainty.mat';
%%  fits every model to each subject data separately 
% cbm_lap(data, @mymodelDDM_twobounds_ORA, prior_Threshold, fname_Threshold);
cbm_lap(data, @mymodel_Optimal_ORA, prior_Optimal, fname_Optimal);
cbm_lap(data, @mymodelUncertainty_ORA, prior_Uncertainty, fname_Uncertainty);
%% check fitted parameter
% fname = load('lap_Threshold.mat');
% cbm   = fname.cbm;
% % look at fitted parameters
% cbm.output.parameters

fname = load('lap_Optimal.mat');
cbm   = fname.cbm;
% look at fitted parameters
cbm.output.parameters

fname = load('lap_Uncertainty.mat');
cbm   = fname.cbm;
% look at fitted parameters
cbm.output.parameters
%% hierarchical Bayesian inference 
% 1st input: data for all subjects
% data
% 2nd input: a cell input containing function handle to models
% models = {@mymodelDDM_twobounds_ORA,@mymodel_Optimal_ORA,@mymodelUncertainty_ORA};
models = {@mymodel_Optimal_ORA,@mymodelUncertainty_ORA};
% 3rd input: another cell input containing file-address to files saved by cbm_lap
% fcbm_maps = {'lap_Threshold.mat','lap_Optimal.mat','lap_Uncertainty.mat'};
fcbm_maps = {'lap_Optimal.mat','lap_Uncertainty.mat'};
% note that they corresponds to models (so pay attention to the order)
% 4th input: a file address for saving the output
fname_hbi = 'hbi_Optimal_Uncertainty.mat';

% run cbm_hbi
cbm_hbi(data,models,fcbm_maps,fname_hbi);
% Running this command, prints a report on your matlab output 
% (e.g. on the command window)
%% take a look at the saved file
fname_hbi = load('hbi_Threhsold_Optimal_Uncertainty.mat');
cbm   = fname_hbi.cbm;
cbm.output

% model frequency
model_frequency = cbm.output.model_frequency

% estimeated group mean
group_mean_Threshold = cbm.output.group_mean{1}
% group mean for parameters of model_RL
group_mean_Optimal = cbm.output.group_mean{2}
% group mean for parameters of model_RLF
group_mean_Uncertainty = cbm.output.group_mean{3}
% group mean for parameters of model_dualRL

% errorbar of estimated parameters
group_errorbar_Threshold = cbm.output.group_hierarchical_errorbar{1};
group_errorbar_Optimal = cbm.output.group_hierarchical_errorbar{2};
group_errorbar_Uncertainty = cbm.output.group_hierarchical_errorbar{3};
%% plot group parameters/main outputs of HBI
% 1st input is the file-address of the file saved by cbm_hbi
fname_hbi = 'hbi_Optimal_Uncertainty.mat';

% 2nd input: a cell input containing model names
model_names = {'Optimal', 'Uncertainty'};
% note that they corresponds to models (so pay attention to the order)

% 3rd input: another cell input containing parameter names of the winning model
param_names = {'\alpha^+','\alpha^+','\alpha^-','\beta'};
% note that '\alpha^+' is in the latex format, which generates a latin alpha

% 4th input: another cell input containing transformation function associated with each parameter of the winning model
transform = {'sigmoid','sigmoid','sigmoid','exp'};
% note that if you use a less usual transformation function, you should pass the handle here (instead of a string)
% cbm_hbi_plot(fname_hbi, model_names, param_names)

cbm_hbi_plot(fname_hbi, model_names, param_names, transform)
% this function creates a model comparison plot (exceednace probability and model frequency) as well as 
% a plot of transformed parameters of the most frequent model.
%% protected exceedance probabilities
close all;clear all
load summaryORA3_TD.mat
data_all = data;
clear data
subjvec = unique(data_all.subjid);
for subjidx = 1:length(subjvec) % for generated data run 20 subjects
    subjidx
    idx = find(data_all.subjid == subjvec(subjidx));
    data{subjidx,1}.red    = data_all.red(idx);
    data{subjidx,1}.green  = data_all.green(idx);
    data{subjidx,1}.choice = data_all.choice(idx);
end

fname_hbi = 'hbi_Optimal_Uncertainty';

% 1st input is data, 
% 2nd input is the file-address of the file saved by cbm_hbi
cbm_hbi_null(data,fname_hbi);
% Running this command, prints a report on your matlab output 
% (e.g. on the command window)

fname_hbi = load('hbi_Optimal_Uncertainty.mat');
cbm   = fname_hbi.cbm;
pxp   = cbm.output.protected_exceedance_prob