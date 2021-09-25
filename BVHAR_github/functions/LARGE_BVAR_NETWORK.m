% LARGE_BVAR_NETWORK

% LINE 33 n=12 THIS ADDS:
% CBOE RUSSELL 2000 VOLATILITY INDEX - PRICE INDEX	
% CBOE DJ INDUSTRIAL VOLATILITY VXD - PRICE INDEX	
% CBOE EMERG. MARKETS VOLATILITY INDEX - PRICE INDEX	
% CBOE NDX VOLATILITY (VXN) - PRICE INDEX	
% CBOE BRAZIL ETF VOLATILITY INDEX - PRICE INDEX	
% CBOE CHINA ETF VOLATILITY INDEX - PRICE INDEX

% to the BVAR-t-CSV-MA(1) model and tracks network connections

% LINE 34 removes adds in the EVZ foreign currency VIX ETF and RUSSELL 2000.

% LINE 35 removes the VXEFA and VIX and replaces with CBOE EMERGING 
% AND CBOE RUSSEL 2000.

% LINE 36 removes the VXEFA and replaces with CBOE BRAZIL and CBOE CHINA
% VIX. It also adds in the EVZ foreign currency VIX ETF.

% LINE 37 USES CBOE VIX INDICES ON:

% GOLDMAN SACHS
% AMAZON
% APPLE
% OVX
% XLE
% GVX
% SLV
% EVZ

clear; clc;
for dataset=1:6
    dataset
    if dataset==1
    data=xlsread('IV_EXT','DATA2','B2:S2473'); 
    data=[data(:,1:11), data(:,14:end)];
    elseif dataset==2
    data=xlsread('IV_EXT','DATA2','B2:I2473'); 
    elseif dataset==3
    data=xlsread('IV_EXT','DATA3','B2:G2473');
    elseif dataset==4
    data=xlsread('IV_EXT','DATA2','B2:K2473'); 
    elseif dataset==5
    data=xlsread('IV_EXT','DATA2','B2:L2473');
    data=[data(:,1:8),data(:,10:end)];
    elseif dataset==6
    data=xlsread('IV_EXT','DAT','B2:T2473');    
    data=[data(:,1:9) data(:,11:13)];
    end
LOG=1;
if LOG==1
data=log(data);
end

window=252; % this is the estimation window that we roll throughout the sample
p=22; nsims=1000; burnin=500;
L=3; % numer of HAR horizons.
cp_ml=0;
model=7; % always model 7 because this is preferred by the data.
IND=(1:1:p+window)';
numsamples=(window+p:1:length(data)); numsamples=length(numsamples); % number of estimations
LIND=length(IND);
% Create some storage matrices
ML_S=single(zeros(numsamples,1));
TF_C=single(zeros(numsamples,nsims));
TC_S=TF_C;
TC_M=TF_C;
TC_L=TF_C;
TN_C=single(zeros(numsamples,size(data,2),nsims));
ND_S=TN_C;
ND_M=TN_C;
ND_L=TN_C;
tic;
%for kk=1:1
for kk=1:numsamples
   data1=data(IND,:);
   Y0=data1(1:p,:);
   shortY=data1(p+1:end,:);
   [T,n]=size(shortY);
   k=n*L+1;
   
   % Storage for connectedness measures
   tfc=zeros(1,nsims);
   tcs=tfc; tcm=tfc; tcl=tfc; 
   ndcl=zeros(n,nsims); ndcm=ndcl; ndcs=ndcl; tndc=ndcl;

   switch model
       case 1
           BVAR_t;
       case 2
           BVAR_CSV;
       case 3
           BVAR_MA;
       case 4
           BVAR_t_CSV;
       case 5
           BVAR_t_MA;
       case 6 
           BVAR_CSV_MA;
       case 7
           BVAR_CSV_t_MA;
   end
   
%   ML_S(kk,1)=ML;
   TC_S(kk,:)=tcs; TC_M(kk,:)=tcm; TC_L(kk,:)=tcl;
   TF_C(kk,:)=tfc; 
   ND_S(kk,:,:)=ndcs; ND_M(kk,:,:)=ndcm; ND_L(kk,:,:)=ndcl;
   TN_C(kk,:,:)=tndc;
   IND=IND+1;  
end
toc; 
%% Save results
switch model
       case 1
           name='LARGE_BVAR_t';
       case 2
           name='LARGE_BVAR_CSV';
       case 3
           name='LARGE_BVAR_MA';
       case 4
           name='LARGE_BVAR_t_CSV';
       case 5
           name='LARGE_BVAR_t_MA';
       case 6 
           name='LARGE_BVAR_CSV_MA';
       case 7
           if LOG==1
           if dataset==1
           name='log_LARGE_BVAR_CSV_t_MA_N18';
           elseif dataset==2
           name='log_LARGE_BVAR_CSV_t_MA_EVZRUT';  
           elseif dataset==3
           name='log_LARGE_BVAR_CSV_t_MA_EEMRUT'; 
           elseif dataset==4
           name='log_LARGE_BVAR_CSV_t_MA_EVZRUTDJEEM';
           elseif dataset==5
           name='log_LARGE_BVAR_CSV_t_MA_EVZRUTEEMND';
           elseif dataset==6
           name='log_LARGE_BVAR_CSV_t_MA_REV';
           end
           else
           if dataset==1
           name='LARGE_BVAR_CSV_t_MA_N18';
           elseif dataset==2
           name='LARGE_BVAR_CSV_t_MA_EVZRUT';  
           elseif dataset==3
           name='LARGE_BVAR_CSV_t_MA_EEMRUT'; 
           elseif dataset==4
           name='LARGE_BVAR_CSV_t_MA_EVZRUTDJEEM';
           elseif dataset==5
           name='LARGE_BVAR_CSV_t_MA_EVZRUTEEMND';
           elseif dataset==6
           name='LARGE_BVAR_CSV_t_MA_REV';
           end
           end
end
%varname(1,:)='ML_S'; varname(2,:)='TC_S'; varname(3,:)='TC_M'; varname(4,:)='TC_L';
%varname(5,:)='TF_C'; varname(6,:)='ND_S'; varname(7,:)='ND_M'; varname(8,:)='ND_L';
%varname(9,:)='TN_C';
varname(1,:)='TC_S'; varname(2,:)='TC_M'; varname(3,:)='TC_L';
varname(4,:)='TF_C'; varname(5,:)='ND_S'; varname(6,:)='ND_M'; varname(7,:)='ND_L';
varname(8,:)='TN_C';
save(name,varname(1,:),varname(2,:),varname(3,:),varname(4,:),varname(5,:),...
    varname(6,:),varname(7,:),varname(8,:));
clearvars
end
