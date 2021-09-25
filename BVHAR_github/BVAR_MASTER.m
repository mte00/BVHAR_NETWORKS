%% Implied Volatility Spillovers of Stocks, Oil and Commodities
%  A Bayesian VAR Analysis with different covariance structures

% Bayesian VARs with Fat Tails using Chan (2020) JBES estimation
% methodology

% Get Frequency connections using Barunik and Krehlik (2018) JFEC

% USE 7 MODELS IN CHAN (2020) ROLLING OVER ANNUAL WINDOW. AT EACH ITERATION
% COMPUTE THE NETWORK MEASURES AND THE MARGINAL LIKELIHOODS.

% On line 90 change model to any integer in 1,2,3,4,5,6,7 to get results
% for each of the 7 BVAR models.

% Code for EJOR submission

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------%
%                                   NOTE                                  %
%-------------------------------------------------------------------------%
% This only provides code to estimate network measures. This is because   %
% the ETF data is propietary and therefore cannot be shared.              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First, we need to compute HAR-type lags ala Corsi (2009) JFEC. By Bubak
% et al. (2011) the VHAR model can be written as VAR(p) with restricted
% parameters. In this case we want lag 1, lag 5 and lag 22 so we have days
% weeks and months.

% Note we can also take log(IV) to get networks in Online Appendix
% Compute networks with standard IV to check for differences. 

clear; clc;
addpath('functions')
addpath('Data')
data=xlsread('ENRGY2','DATA','B2:G2473');
LOG=0;
if LOG==1
data=log(data); % comment out if want to use raw data
end

% Get descriptive stats;
STATS=zeros(5,size(data,2));
for kk=1:5
if kk==1    
temp1=mean(data);
elseif kk==2
temp1=median(data);
elseif kk==3
temp1=std(data,0,1);
elseif kk==4
temp1=skewness(data);
elseif kk==5
temp1=kurtosis(data);
end
STATS(kk,:)=temp1; % This is Table 1 of paper
end
CORR=corr(data);
TAB1 = [STATS; CORR]; % This is table 1 of paper
%--------------------------------------------------------------------------
% assess significance of correlation at 1% level using Cryer and Chan 2008
% adjustment.

% First estimate AR(1) models of IV indices to get AR(1) coefficients. 
yy=data(2:end,:); xx=data(1:end-1,:);
for kk=1:size(yy,2)
y=yy(:,kk);
x=[ones(length(xx),1) xx(:,kk)]; 
temp=(x'*x)\x'*y;
b(:,kk)=temp;
end
b=b(2,:);
nn=sqrt(length(yy));
scale=eye(kk);
for kk=1:size(yy,2)
for jj=1:size(yy,2)
    if kk ~=jj
    scale(kk,jj)=(2.58/nn)*sqrt((1+b(kk)*b(jj))/(1-b(kk)*b(jj)));
    end
end
end
ee=CORR>scale; % confirms off diagonals in CORR are greater than scale that contains the 1% value.
clear temp1 yy xx b nn
%--------------------------------------------------------------------------
% BVHAR Estimation 
window=252; % this is the estimation window that we roll throughout the sample
p=22; nsims=1000; burnin=500;
L=3; % number of HAR horizon lags
cp_ml=1;
model=1 % change this to estimate all seven different models. 
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
%--------------------------------------------------------------------------
estimate=0; % need to estimate the models to get results, set this to 1.
if estimate==1
tic;
%for kk=1:1 <#OK#>
for kk=1:numsamples
   data1=data(IND,:);
   Y0=data1(1:p,:);
   shortY=data1(p+1:end,:);
   [T,n]=size(shortY);
   % Since we use HAR-type lags here and only have 3 "horizons" k=n*3+1;
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
           BVAR_CSV_t_MA; % these scripts contain the scaling adjustments to the 
           % covariance matrices as in the paper.
   end
   
   ML_S(kk,1)=ML;
   TC_S(kk,:)=tcs; TC_M(kk,:)=tcm; TC_L(kk,:)=tcl;
   TF_C(kk,:)=tfc; 
   ND_S(kk,:,:)=ndcs; ND_M(kk,:,:)=ndcm; ND_L(kk,:,:)=ndcl;
   TN_C(kk,:,:)=tndc;
   IND=IND+1;  
end
toc; 
% Save results
if LOG==1
switch model
       case 1
           name='log_BVAR_t';
       case 2
           name='log_BVAR_CSV';
       case 3
           name='log_BVAR_MA';
       case 4
           name='log_BVAR_t_CSV';
       case 5
           name='log_BVAR_t_MA';
       case 6 
           name='log_BVAR_CSV_MA';
       case 7
           name='log_BVAR_CSV_t_MA';
end
else
    switch model
       case 1
           name='BVAR_t';
       case 2
           name='BVAR_CSV';
       case 3
           name='BVAR_MA';
       case 4
           name='BVAR_t_CSV';
       case 5
           name='BVAR_t_MA';
       case 6 
           name='BVAR_CSV_MA';
       case 7
           name='BVAR_CSV_t_MA';
    end
end
varname(1,:)='ML_S'; varname(2,:)='TC_S'; varname(3,:)='TC_M'; varname(4,:)='TC_L';
varname(5,:)='TF_C'; varname(6,:)='ND_S'; varname(7,:)='ND_M'; varname(8,:)='ND_L';
varname(9,:)='TN_C';
save(name,varname(1,:),varname(2,:),varname(3,:),varname(4,:),varname(5,:),...
    varname(6,:),varname(7,:),varname(8,:),varname(9,:));
    
elseif estimate==0
end

%% Now to get Results
% Note below WILL NOT WORK unless you have estimated all 7 variants of the
% BVAR. You should also run the below sections separately since some clear
% the workspace.
BVAR_model_selection; % This will give you Table 2 and Figures 1 and A1 of the paper
%%
Network_Analysis; % This will provide the network results and generate 
                  % Figures 2, 3, 4, 5 and A.2.  Note this 
                  % needs to be ran for LOG=1 and LOG=0 to get all
                  % results in paper. When LOG=1 you get corresponding figures in
                  % Appendix
LARGE_BVAR_NETWORK; % Gets Network measures for other datasets. Note this 
                    % needs to be ran for LOG=1 and LOG=0 to get all
                    % results in paper. 
LARGE_Network_Analysis; % LOG=0 and LOG=1 gives you results in Appendix    
%--------------------------------------------------------------------------