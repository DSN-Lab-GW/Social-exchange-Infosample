load data_summary_ASD.mat %_TD/_ASD_matched/_TD_mathed

data_sum{1,i_s}.data(:,1)  %id
data_sum{1,i_s}.data(:,2)  %trial
data_sum{1,i_s}.data(:,3)  %prob
data_sum{1,i_s}.data(:,4)  %choice
data_sum{1,i_s}.data(:,5)  %green
data_sum{1,i_s}.data(:,6)  %red
data_sum{1,i_s}.data(:,7)  %covered tiles
data_sum{1,i_s}.data(:,9)  %0 for red
data_sum{1,i_s}.data(:,9)  %1 for green
data_sum{1,i_s}.data(:,8)  %recoded choice
%% reciprocation per probability
ids = [1:size(data_sum,2)];
prob_s = unique(data_sum{1,1}.data(:,3));
for   i_s =1:length(ids)
    for i_p = 1:length(prob_s)

        prob = find(data_sum{1,i_s}.data(:,3)==prob_s(i_p));
        mean_prob(i_s,i_p) = nanmean(data_sum{1,i_s}.data([prob'],8))
    end 
end

for   i_s =1:length(ids)
    for i_p = 1:length(prob_s)
      
        prob = find(data_sum{1,i_s}.data(:,3)==prob_s(i_p));
        mean_prob_all(10*(i_s-1)+1:10*i_s,i_p) = data_sum{1,i_s}.data([prob'],8);
    end 
end
%% tiles coverd per probability
prob_s = unique(data_sum{1,1}.data(:,3));
for   i_s =1:length(ids)
for i_p = 1:length(prob_s)
 
    tiles_num = find(data_sum{1,i_s}.data(:,3)==prob_s(i_p));
    mean_tiles(i_s,i_p) = nanmean(data_sum{1,i_s}.data([tiles_num'],7))
end 
end

for   i_s =1:length(ids)
    for i_p = 1:length(prob_s)

        tiles_num = find(data_sum{1,i_s}.data(:,3)==prob_s(i_p));
        mean_tiles_all(10*(i_s-1)+1:10*i_s,i_p) = data_sum{1,i_s}.data([tiles_num'],7);
    end 
end

%% prepare for lme in R
for   i_s =1:length(ids)
    for i_p = 1:length(prob_s)
        
        prob = find(data_sum{1,i_s}.data(:,3)==prob_s(i_p));
        mean_diff_all(10*(i_s-1)+1:10*i_s,i_p) = abs(data_sum{1,i_s}.data([prob'],6) - data_sum{1,i_s}.data([prob'],5));
        subjnum4lme(10*(i_s-1)+1:10*i_s,i_p) = i_s;
        age4lme(10*(i_s-1)+1:10*i_s,i_p) = demographic.age(i_s);
        recipro4lme(10*(i_s-1)+1:10*i_s,i_p) = prob_s(i_p);
        srs4lme(10*(i_s-1)+1:10*i_s,i_p) = demographic.srs_raw(i_s);
%         scq4lme(10*(i_s-1)+1:10*i_s,i_p) = demographic.scq(i_s);
%         aq4lme(10*(i_s-1)+1:10*i_s,i_p) = demographic.aq(i_s);
        gender4lme(10*(i_s-1)+1:10*i_s,i_p) = demographic.gender(i_s);
    end 
end
lmedata(:,1) = reshape(subjnum4lme',length(ids)*60,1);% subnum*60
lmedata(:,2) = reshape(mean_tiles_all',length(ids)*60,1);
lmedata(:,3) = reshape(recipro4lme',length(ids)*60,1);
lmedata(:,4) = reshape(mean_prob_all',length(ids)*60,1);
lmedata(:,5) = reshape(mean_diff_all',length(ids)*60,1);
lmedata(:,6) = reshape(age4lme',length(ids)*60,1);
lmedata(:,7) = reshape(srs4lme',length(ids)*60,1);
lmedata(:,8) = reshape(scq4lme',length(ids)*60,1);
lmedata(:,9) = reshape(aq4lme',length(ids)*60,1);
lmedata(:,10) = reshape(gender4lme',length(ids)*60,1);
%% age group plot
age2 = reshape(repmat(demographic.age,1,10)',480,1);%number of subject * 10
A=mink(age2,160) %120
B=maxk(age2,160) %120
agegroup{1} = find(age2<=max(A));
agegroup{2} = find(max(A)<age2 & age2<=min(B));
agegroup{3} = find(age2>min(B));
% invest prob
figure(1)
for i = 1:size(agegroup,2)
    plot(nanmean(mean_prob_all(agegroup{i},:),1))
    xlabel('Probability of reciprocation')
    hold on
    grid off
    xlim([0 7]);
    xticks([0 1 2 3 4 5 6 7 8])
    xticklabels({' ','0','0.2','0.4','0.6','0.8','1',' '})
    ylim([0 1]);
    ylabel('Proportion of invest')

end
title(legend,'Age');
legend('8~9','9~11','11~12','Location','southeast')
for i = 1:size(agegroup,2)
    er = errorbar(nanmean(mean_prob_all(agegroup{i},:),1), (nanstd(mean_prob_all(agegroup{i},:)))/sqrt(length(mean_prob_all(agegroup{i},:))));
    er.Color = [0 0 0];
    er.LineStyle = 'none'; 
    er.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
hold off
% info sample
figure(2)
for i = 1:size(agegroup,2)
    plot(nanmean(mean_tiles_all(agegroup{i},:),1))
    xlabel('Probability of reciprocation')
    hold on
    grid off
    xlim([0 7]);
    xticks([0 1 2 3 4 5 6 7 8])
    xticklabels({' ','0','0.2','0.4','0.6','0.8','1',' '})
    ylim([6 22]);
    ylabel('Num of tiles sampled')

end
title(legend,'Age');
legend('8~9','9~11','11~12','Location','southeast')
for i = 1:size(agegroup,2)
    er = errorbar(nanmean(mean_tiles_all(agegroup{i},:),1), (nanstd(mean_tiles_all(agegroup{i},:)))/sqrt(length(mean_tiles_all(agegroup{i},:))));
    er.Color = [0 0 0];
    er.LineStyle = 'none'; 
    er.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
hold off

%% plot age\sample frequency TD
% plot age\sample frequency (age)
    figure(5)
    scatter(demographic.age/12,nanmean(mean_tiles,2))
    ylim([0 25])
    ylabel('Num of tiles sampled')
    xlim([7 13])
    xlabel('Age(month)')
    hold on
    Fit = polyfit(demographic.age/12,nanmean(mean_tiles,2),2); % x = x data, y = y data, 1 = order of the polynomial i.e a straight line 
    % Make 100 fitted samples going from -13 to +12.
    xFit = linspace(1, 20, 100);
    % Get the estimated values with polyval()
    yFit = polyval(Fit, xFit);
    % Plot the fit
    plot(xFit, yFit, 'b-', 'LineWidth', 2);
    hold off
    
%% TD&ASD tiles plot
load data_summary_TD_20220403.mat mean_tiles_all
% load data_summary_TD_male.mat mean_tiles_all
mean_tilesTD = mean_tiles_all;
clear mean_tiles_all
load data_summary_ASD_20220831.mat mean_tiles_all
mean_tilesASD = mean_tiles_all;
clear mean_tiles_all

figure(3)
plot(nanmean(mean_tilesTD,1))%,'color','[0.8500 0.3250 0.0980]'
hold on
grid off

plot(nanmean(mean_tilesASD,1))%,'color','[0.8500 0.3250 0.0980]'
xlabel('Probability of reciprocation')
xlim([0 7]);
xticks([0 1 2 3 4 5 6 7 8])
xticklabels({' ','0','0.2','0.4','0.6','0.8','1',' '})
ylim([0 25]);
ylabel('Num of tiles sampled')
legend('TD','ASD','Location','southeast')
er = errorbar(nanmean(mean_tilesTD,1), (nanstd(mean_tilesTD))/sqrt(length(mean_tilesTD)));
er.Color = [0 0 0];
er.LineStyle = 'none'; 
er.Annotation.LegendInformation.IconDisplayStyle = 'off';

er = errorbar(nanmean(mean_tilesASD,1), (nanstd(mean_tilesASD))/sqrt(length(mean_tilesASD)));
er.Color = [0 0 0];
er.LineStyle = 'none'; 
er.Annotation.LegendInformation.IconDisplayStyle = 'off';
hold off
%% TD&ASD invest plot
load data_summary_TD_male.mat mean_prob_all
mean_probTD = mean_prob_all;
clear mean_prob_all
load data_summary_ASD_model_male.mat mean_prob_all
mean_probASD = mean_prob_all;
clear mean_prob_all

figure(4)
plot(nanmean(mean_probTD,1))%,'color','[0.8500 0.3250 0.0980]'
hold on
grid off

plot(nanmean(mean_probASD,1))%,'color','[0.8500 0.3250 0.0980]'
xlabel('Probability of reciprocation')
xlim([0 7]);
xticks([0 1 2 3 4 5 6 7 8])
xticklabels({' ','0','0.2','0.4','0.6','0.8','1',' '})
ylim([0 1]);
ylabel('Proportion of invest')
legend('TD','ASD','Location','southeast')
er = errorbar(nanmean(mean_probTD,1), (nanstd(mean_probTD))/sqrt(length(mean_probTD)));
er.Color = [0 0 0];
er.LineStyle = 'none'; 
er.Annotation.LegendInformation.IconDisplayStyle = 'off';

er = errorbar(nanmean(mean_probASD,1), (nanstd(mean_probASD))/sqrt(length(mean_probASD)));
er.Color = [0 0 0];
er.LineStyle = 'none'; 
er.Annotation.LegendInformation.IconDisplayStyle = 'off';
hold off


%% demographic diff test
[h,p,ci,stats] = ttest2(demographicASD.age,demographicTD.age)
[h,p,ci,stats] = ttest2(demographicASD.srs_raw,demographicTD.srs_raw)
[h,p,ci,stats] = ttest2(demographicASD.age,demographicTD.aq)
[h,p,ci,stats] = ttest2(demographicASD.age,demographicTD.scq)
hist(demographicASD.scq)
hist(demographicASD.srs_raw)
V = nanvar(demographicASD.scq)
V = nanvar(demographicASD.srs_raw)
[R,P] = corr(demographicTD.aq,demographicTD.srs_raw,'rows','complete')

err = [ nanstd(demographicTD.age)/sqrt(length(demographicTD.age)) nanstd(demographicASD.age)/sqrt(length(demographicASD.age));
        nanstd(demographicTD.srs_raw) nanstd(demographicASD.srs_raw);
        nanstd(demographicTD.aq) nanstd(demographicASD.aq);
        nanstd(demographicTD.scq) nanstd(demographicASD.scq)]
%chi squre
 % Observed data
       n1 = sum(demographicTD.gender); N1 = size(demographicTD.gender,1);
       n2 = sum(demographicASD.gender); N2 = size(demographicASD.gender,1);
       % Pooled estimate of proportion
       p0 = (n1+n2) / (N1+N2)
       % Expected counts under H0 (null hypothesis)
       n10 = N1 * p0;
       n20 = N2 * p0;
       % Chi-square test, by hand
       observed = [n1 N1-n1 n2 N2-n2];
       expected = [n10 N1-n10 n20 N2-n20];
       chi2stat = sum((observed-expected).^2 ./ expected)
       p = 1 - chi2cdf(chi2stat,1)
%% exclued subject sample all the time/no sample at all
load data_summary_TD_20220403.mat
load data_summary_ASD_20220420.mat

morethan90lessthan10 = find(sum(mean_tiles,2)<=25*6*0.10 | sum(mean_tiles,2)>=25*6*0.90)
for i = 1:length(morethan90lessthan10)
    data(data(:,1) == morethan90lessthan10(i),:) = [];
    sub_code{morethan90lessthan10(i)} = [];
    data_sum{morethan90lessthan10(i)} = [];
    demographic.age(morethan90lessthan10(i)) = 0;
    demographic.srs_raw(morethan90lessthan10(i)) = 0;
    demographic.aq(morethan90lessthan10(i)) = 0;
    demographic.scq(morethan90lessthan10(i)) = 1000;
    demographic.gender(morethan90lessthan10(i)) = nan;
end

sub_code = sub_code(~cellfun('isempty',sub_code));
data_sum = data_sum(~cellfun('isempty',data_sum));
demographicTD.age = nonzeros(demographic.age);
demographicTD.srs_raw = nonzeros(demographic.srs_raw);
demographicTD.scq = demographic.scq(demographic.scq~=1000);
demographicTD.aq = nonzeros(demographic.aq);
demographicTD.gender = (demographic.gender(~isnan(demographic.gender)));
data_all = data;
clear data;clear demographic

sub_code = sub_code(~cellfun('isempty',sub_code));
data_sum = data_sum(~cellfun('isempty',data_sum));
demographicASD.age = nonzeros(demographic.age);
demographicASD.srs_raw = nonzeros(demographic.srs_raw);
demographicASD.scq = demographic.scq(demographic.scq~=1000);
demographicASD.aq = nonzeros(demographic.aq);
demographicASD.gender = (demographic.gender(~isnan(demographic.gender)));
data_all = data;
clear data;clear demographic
%% make model format
data.subjid = data_all(:,1); %id
data.trial = data_all(:,3); %trial
data.recip = data_all(:,4); %prob
data.rawChoice = data_all(:,5); %choice
data.green = data_all(:,6); %green
data.red = data_all(:,7); %red
data.closed = data_all(:,8); %covered tiles
data.choice(data.rawChoice == 0,1) = 1;
data.choice(data.rawChoice ~= 0,1) = -1;
data.open = 25-data.closed;
demographic = demographicTD;
% demographic = demographicASD;



