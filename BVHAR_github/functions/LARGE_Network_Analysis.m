% Network Analysis LARGE BVAR
clear;clc;
addpath('Data'); addpath('functions')
% Import dates
opts = delimitedTextImportOptions("NumVariables", 2);
LOG=0;
% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["VarName1", "x"];
opts.VariableTypes = ["double", "datetime"];
opts = setvaropts(opts, 2, "InputFormat", "dd/MM/yyyy");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

Dates = readtable("Dates.csv", opts);
dates2=Dates.x;
dates2=dates2(2:end,:);
dates1=Dates.x;
dates1=dates1(252+22+1:end);
datetime.setDefaultFormats('defaultdate','yyyy-MM-dd')
datnum=datenum(dates1);
dates=datetime(datnum,'ConvertFrom','datenum','Format','yyyy-MM-dd');

numsamples=2199;
% Horizon Specific Network Connectedness
temp = zeros(numsamples,5);
temp2 = zeros(numsamples,5);
temp3 = zeros(numsamples,5);
% load
if LOG==0
load('LARGE_BVAR_CSV_t_MA_N18', 'TC_S')
temp(:,1)=median(TC_S,2);
% load
load('LARGE_BVAR_CSV_t_MA_EVZRUT', 'TC_S')
temp(:,2)=median(TC_S,2);
% load
load('LARGE_BVAR_CSV_t_MA_EEMRUT', 'TC_S')
temp(:,3)=median(TC_S,2);
% load
load('LARGE_BVAR_CSV_t_MA_EVZRUTDJEEM', 'TC_S')
temp(:,4)=median(TC_S,2);
% load
load('LARGE_BVAR_CSV_t_MA_REV', 'TC_S')
temp(:,5)=median(TC_S,2);
% load
load('LARGE_BVAR_CSV_t_MA_N18', 'TC_M')
temp2(:,1)=median(TC_M,2);
% load
load('LARGE_BVAR_CSV_t_MA_EVZRUT', 'TC_M')
temp2(:,2)=median(TC_M,2);
% load
load('LARGE_BVAR_CSV_t_MA_EEMRUT', 'TC_M')
temp2(:,3)=median(TC_M,2);
% load
load('LARGE_BVAR_CSV_t_MA_EVZRUTDJEEM', 'TC_M')
temp2(:,4)=median(TC_M,2);
% load
load('LARGE_BVAR_CSV_t_MA_REV', 'TC_M')
temp2(:,5)=median(TC_M,2);
% load
load('LARGE_BVAR_CSV_t_MA_N18', 'TC_L')
temp3(:,1)=median(TC_L,2);
% load
load('LARGE_BVAR_CSV_t_MA_EVZRUT', 'TC_L')
temp3(:,2)=median(TC_L,2);
% load
load('LARGE_BVAR_CSV_t_MA_EEMRUT', 'TC_L')
temp3(:,3)=median(TC_L,2);
% load
load('LARGE_BVAR_CSV_t_MA_EVZRUTDJEEM', 'TC_L')
temp3(:,4)=median(TC_L,2);
% load
load('LARGE_BVAR_CSV_t_MA_REV', 'TC_L')
temp3(:,5)=median(TC_L,2);
elseif LOG==1
load('log_LARGE_BVAR_CSV_t_MA_N18', 'TC_S')
temp(:,1)=median(TC_S,2);
% load
load('log_LARGE_BVAR_CSV_t_MA_EVZRUT', 'TC_S')
temp(:,2)=median(TC_S,2);
% load
load('log_LARGE_BVAR_CSV_t_MA_EEMRUT', 'TC_S')
temp(:,3)=median(TC_S,2);
% load
load('log_LARGE_BVAR_CSV_t_MA_EVZRUTDJEEM', 'TC_S')
temp(:,4)=median(TC_S,2);
% load
load('log_LARGE_BVAR_CSV_t_MA_REV', 'TC_S')
temp(:,5)=median(TC_S,2);
% load
load('log_LARGE_BVAR_CSV_t_MA_N18', 'TC_M')
temp2(:,1)=median(TC_M,2);
% load
load('log_LARGE_BVAR_CSV_t_MA_EVZRUT', 'TC_M')
temp2(:,2)=median(TC_M,2);
% load
load('log_LARGE_BVAR_CSV_t_MA_EEMRUT', 'TC_M')
temp2(:,3)=median(TC_M,2);
% load
load('log_LARGE_BVAR_CSV_t_MA_EVZRUTDJEEM', 'TC_M')
temp2(:,4)=median(TC_M,2);
% load
load('log_LARGE_BVAR_CSV_t_MA_REV', 'TC_M')
temp2(:,5)=median(TC_M,2);
% load
load('log_LARGE_BVAR_CSV_t_MA_N18', 'TC_L')
temp3(:,1)=median(TC_L,2);
% load
load('log_LARGE_BVAR_CSV_t_MA_EVZRUT', 'TC_L')
temp3(:,2)=median(TC_L,2);
% load
load('log_LARGE_BVAR_CSV_t_MA_EEMRUT', 'TC_L')
temp3(:,3)=median(TC_L,2);
% load
load('log_LARGE_BVAR_CSV_t_MA_EVZRUTDJEEM', 'TC_L')
temp3(:,4)=median(TC_L,2);
% load
load('log_LARGE_BVAR_CSV_t_MA_REV', 'TC_L')
temp3(:,5)=median(TC_L,2);
end


X=table(dates,temp(:,1),temp(:,2),temp(:,3),temp(:,4),temp(:,5));
if LOG==0
writetable(X,'Short2.csv')
elseif LOG==1
writetable(X,'Short2_log.csv')   
end
X=table(dates,temp2(:,1),temp2(:,2),temp2(:,3),temp2(:,4),temp2(:,5));
if LOG==0
writetable(X,'Medium2.csv')
elseif LOG==1
writetable(X,'Medium2_log.csv')   
end
X=table(dates,temp3(:,1),temp3(:,2),temp3(:,3),temp3(:,4),temp3(:,5));
if LOG==0
writetable(X,'Long2.csv')
elseif LOG==1
writetable(X,'Long2_log.csv')   
end

figure(22)
subplot(3,1,1)
plot(temp(:,1),'b','LineWidth',1.5)
hold on,
plot(temp(:,2),'r','LineWidth',1.5)
plot(temp(:,3),'m','LineWidth',1.5)
plot(temp(:,4),'k','LineWidth',1.5)
plot(temp(:,5),'g','LineWidth',1.5)
axis tight
grid on
ylabel('Short-term')
subplot(3,1,2)
plot(temp2(:,1),'b','LineWidth',1.5)
hold on,
plot(temp2(:,2),'r','LineWidth',1.5)
plot(temp2(:,3),'m','LineWidth',1.5)
plot(temp2(:,4),'k','LineWidth',1.5)
plot(temp2(:,5),'g','LineWidth',1.5)
axis tight
grid on
ylabel('Medium-term')
subplot(3,1,3)
plot(temp3(:,1),'b','LineWidth',1.5)
hold on,
plot(temp3(:,2),'r','LineWidth',1.5)
plot(temp3(:,3),'m','LineWidth',1.5)
plot(temp3(:,4),'k','LineWidth',1.5)
plot(temp3(:,5),'g','LineWidth',1.5)
ylabel('Long-term')
axis tight
grid on
legend('Model 1','Model 2', 'Model 3',...
    'Model 4','Model 5','Location','NorthWest')



