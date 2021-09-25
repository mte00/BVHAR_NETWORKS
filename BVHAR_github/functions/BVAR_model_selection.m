% GET BAYES FACTORS
%clear; clc;
data=xlsread('ENRGY2','DATA','B2:G2473');
window=252; % this is the estimation window that we roll throughout the sample
p=22; nsims=1000;
numsamples=(window+p:1:length(data)); numsamples=length(numsamples); % number of estimations
temp = zeros(numsamples,7);

if LOG==1
% load in BVAR-t
load('log_BVAR_t.mat', 'ML_S')
temp(:,1)=ML_S;
% load in BVAR-CSV
load('log_BVAR_CSV.mat', 'ML_S')
temp(:,2)=ML_S; % this has weird values toward end of sample
% load in BVAR-MA
load('log_BVAR_MA.mat', 'ML_S')
temp(:,3)=ML_S;
% load in BVAR-CSV-t
load('log_BVAR_t_CSV.mat', 'ML_S')
temp(:,4)=ML_S;
% load in BVAR-t-MA
load('log_BVAR_t_MA.mat', 'ML_S')
temp(:,5)=ML_S;
% load in BVAR-CSV-MA
load('log_BVAR_CSV_MA.mat', 'ML_S')
temp(:,6)=ML_S; % this has weird values toward end of sample
% load in BVAR-CSV-t-MA
load('log_BVAR_CSV_t_MA.mat', 'ML_S')
temp(:,7)=ML_S;
clear ML_S   
else
% load in BVAR-t
load('BVAR_t.mat', 'ML_S')
temp(:,1)=ML_S;
% load in BVAR-CSV
load('BVAR_CSV.mat', 'ML_S')
temp(:,2)=ML_S; % this has weird values toward end of sample
% load in BVAR-MA
load('BVAR_MA.mat', 'ML_S')
temp(:,3)=ML_S;
% load in BVAR-CSV-t
load('BVAR_t_CSV.mat', 'ML_S')
temp(:,4)=ML_S;
% load in BVAR-t-MA
load('BVAR_t_MA.mat', 'ML_S')
temp(:,5)=ML_S;
% load in BVAR-CSV-MA
load('BVAR_CSV_MA.mat', 'ML_S')
temp(:,6)=ML_S; % this has weird values toward end of sample
% load in BVAR-CSV-t-MA
load('BVAR_CSV_t_MA.mat', 'ML_S')
temp(:,7)=ML_S;
clear ML_S
end
%
% BVAR-CSV and BVAR-CSV-MA(1) were models that gave warnings around badly
% scaled matrices for the log marginal likelihoods toward end of sample.
% Interestingly both have normal errors along with the common stochastic
% volatility term. Get some weird results using log RVs.... however LBFs
% indicate to choose the BVAR_CSV_t_MA models

% We remove these outliers by inserting NaN and then filling the missing
% data with the 1-month moving median value.

% threshold for BVAR-CSV is -3316.65 which is obs 2216, day T
% threshold for BVAR-CSV_MA is -3293.77 which is obs 2215, day T-1

A=temp(:,2);
%A(A<temp(end,2))=NaN;
B=temp(:,6);
%B(B<temp(end-1,6))=NaN;
%A=fillmissing(A,'movmedian',20);
%B=fillmissing(B,'movmedian',20);

temp=[temp(:,1) A, temp(:,3:5), B, temp(:,end)];
clear A B

AVE_ML = mean(temp);
MED_ML = median(temp);

TAB2=[AVE_ML; MED_ML]'; % This is Table 2 of paper

% BVAR-t-CSV-MA(1) has the highest average and median ML so now compute
% Bayes factors and plot for all models!

temp2=cumsum(temp);

temp3=zeros(length(temp),7-1);
for kk=1:6
   temp3(:,kk)=temp2(:,end)-temp2(:,kk);
end

figure(1)
plot(temp3(:,1),'b','LineWidth',1.5)
hold on,
plot(temp3(:,2),'r','LineWidth',1.5)
plot(temp3(:,3),'m','LineWidth',1.5)
plot(temp3(:,4),'k','LineWidth',1.5)
plot(temp3(:,5),'g','LineWidth',1.5)
plot(temp3(:,6),'c','LineWidth',1.5)
legend('BVAR-t-CSV-MA(1) vs. BVAR-t','BVAR-t-CSV-MA(1) vs. BVAR-CSV', 'BVAR-t-CSV-MA(1) vs. BVAR-MA(1)',...
    'BVAR-t-CSV-MA(1) vs. BVAR-t-CSV','BVAR-t-CSV-MA(1) vs. BVAR-t-MA(1)','BVAR-t-CSV-MA(1) vs. BVAR-CSV-MA(1)','Location','NorthWest')
axis tight
grid on
%%
% Import dates
opts = delimitedTextImportOptions("NumVariables", 2);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["VarName1", "x"];
opts.VariableTypes = ["double", "datetime"];
opts = setvaropts(opts, 2, "InputFormat", "dd/MM/yyyy");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
Dates = readtable("Dates.csv", opts);
dates2=Dates.x;
dates2=dates2(2:end,:);
dates1=Dates.x;
dates1=dates1(252+p+1:end);
datetime.setDefaultFormats('defaultdate','yyyy-MM-dd')
datnum=datenum(dates1);
dates=datetime(datnum,'ConvertFrom','datenum','Format','yyyy-MM-dd');
DD=dates; % hold dates here for netowrk plots.
datnum2=datenum(dates2);
dates2=datetime(datnum2,'ConvertFrom','datenum','Format','yyyy-MM-dd');

temp=table(dates,temp3(:,1),temp3(:,2),temp3(:,3),temp3(:,4),temp3(:,5),temp3(:,6));
writetable(temp,'LBF.csv')
% Writes to csv files so can plot as tikzpicture

% Plot Data for Appendix
dates=dates2;
temp=table(dates,data(:,1),data(:,2),data(:,3),data(:,4),data(:,5),data(:,6));
writetable(temp,'DATA.csv')

dates=DD;



