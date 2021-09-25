% Network Analysis
numsamples=2199;
% Horizon Specific Network Connectedness
temp = zeros(numsamples,7);
temp2 = zeros(numsamples,7);
temp3 = zeros(numsamples,7);


if LOG==0
% load in BVAR-t
load('BVAR_t.mat', 'TC_S')
temp(:,1)=median(TC_S,2);
% load in BVAR-CSV
load('BVAR_CSV.mat', 'TC_S')
temp(:,2)=median(TC_S,2); % this has weird values toward end of sample
% load in BVAR-MA
load('BVAR_MA.mat', 'TC_S')
temp(:,3)=median(TC_S,2);
% load in BVAR-CSV-t
load('BVAR_t_CSV.mat', 'TC_S')
temp(:,4)=median(TC_S,2);
% load in BVAR-t-MA
load('BVAR_t_MA.mat', 'TC_S')
temp(:,5)=median(TC_S,2);
% load in BVAR-CSV-MA
load('BVAR_CSV_MA.mat', 'TC_S')
temp(:,6)=median(TC_S,2); % this has weird values toward end of sample
% load in BVAR-CSV-t-MA
load('BVAR_CSV_t_MA.mat', 'TC_S')
temp(:,7)=median(TC_S,2);
clear TC_S

% load in BVAR-t
load('BVAR_t.mat', 'TC_M')
temp2(:,1)=median(TC_M,2);
% load in BVAR-CSV
load('BVAR_CSV.mat', 'TC_M')
temp2(:,2)=median(TC_M,2); % this has weird values toward end of sample
% load in BVAR-MA
load('BVAR_MA.mat', 'TC_M')
temp2(:,3)=median(TC_M,2);
% load in BVAR-CSV-t
load('BVAR_t_CSV.mat', 'TC_M')
temp2(:,4)=median(TC_M,2);
% load in BVAR-t-MA
load('BVAR_t_MA.mat', 'TC_M')
temp2(:,5)=median(TC_M,2);
% load in BVAR-CSV-MA
load('BVAR_CSV_MA.mat', 'TC_M')
temp2(:,6)=median(TC_M,2); % this has weird values toward end of sample
% load in BVAR-CSV-t-MA
load('BVAR_CSV_t_MA.mat', 'TC_M')
temp2(:,7)=median(TC_M,2);
clear TC_M

% load in BVAR-t
load('BVAR_t.mat', 'TC_L')
temp3(:,1)=median(TC_L,2);
% load in BVAR-CSV
load('BVAR_CSV.mat', 'TC_L')
temp3(:,2)=median(TC_L,2); % this has weird values toward end of sample
% load in BVAR-MA
load('BVAR_MA.mat', 'TC_L')
temp3(:,3)=median(TC_L,2);
% load in BVAR-CSV-t
load('BVAR_t_CSV.mat', 'TC_L')
temp3(:,4)=median(TC_L,2);
% load in BVAR-t-MA
load('BVAR_t_MA.mat', 'TC_L')
temp3(:,5)=median(TC_L,2);
% load in BVAR-CSV-MA
load('BVAR_CSV_MA.mat', 'TC_L')
temp3(:,6)=median(TC_L,2); % this has weird values toward end of sample
% load in BVAR-CSV-t-MA
load('BVAR_CSV_t_MA.mat', 'TC_L')
temp3(:,7)=median(TC_L,2);
clear TC_L
elseif LOG==1
% load in log_BVAR-t
load('log_BVAR_t.mat', 'TC_S')
temp(:,1)=median(TC_S,2);
% load in log_BVAR-CSV
load('log_BVAR_CSV.mat', 'TC_S')
temp(:,2)=median(TC_S,2); % this has weird values toward end of sample
% load in log_BVAR-MA
load('log_BVAR_MA.mat', 'TC_S')
temp(:,3)=median(TC_S,2);
% load in log_BVAR-CSV-t
load('log_BVAR_t_CSV.mat', 'TC_S')
temp(:,4)=median(TC_S,2);
% load in log_BVAR-t-MA
load('log_BVAR_t_MA.mat', 'TC_S')
temp(:,5)=median(TC_S,2);
% load in log_BVAR-CSV-MA
load('log_BVAR_CSV_MA.mat', 'TC_S')
temp(:,6)=median(TC_S,2); % this has weird values toward end of sample
% load in log_BVAR-CSV-t-MA
load('log_BVAR_CSV_t_MA.mat', 'TC_S')
temp(:,7)=median(TC_S,2);
clear TC_S

% load in log_BVAR-t
load('log_BVAR_t.mat', 'TC_M')
temp2(:,1)=median(TC_M,2);
% load in log_BVAR-CSV
load('log_BVAR_CSV.mat', 'TC_M')
temp2(:,2)=median(TC_M,2); % this has weird values toward end of sample
% load in log_BVAR-MA
load('log_BVAR_MA.mat', 'TC_M')
temp2(:,3)=median(TC_M,2);
% load in log_BVAR-CSV-t
load('log_BVAR_t_CSV.mat', 'TC_M')
temp2(:,4)=median(TC_M,2);
% load in log_BVAR-t-MA
load('log_BVAR_t_MA.mat', 'TC_M')
temp2(:,5)=median(TC_M,2);
% load in log_BVAR-CSV-MA
load('log_BVAR_CSV_MA.mat', 'TC_M')
temp2(:,6)=median(TC_M,2); % this has weird values toward end of sample
% load in log_BVAR-CSV-t-MA
load('log_BVAR_CSV_t_MA.mat', 'TC_M')
temp2(:,7)=median(TC_M,2);
clear TC_M

% load in log_BVAR-t
load('log_BVAR_t.mat', 'TC_L')
temp3(:,1)=median(TC_L,2);
% load in log_BVAR-CSV
load('log_BVAR_CSV.mat', 'TC_L')
temp3(:,2)=median(TC_L,2); % this has weird values toward end of sample
% load in log_BVAR-MA
load('log_BVAR_MA.mat', 'TC_L')
temp3(:,3)=median(TC_L,2);
% load in log_BVAR-CSV-t
load('log_BVAR_t_CSV.mat', 'TC_L')
temp3(:,4)=median(TC_L,2);
% load in log_BVAR-t-MA
load('log_BVAR_t_MA.mat', 'TC_L')
temp3(:,5)=median(TC_L,2);
% load in log_BVAR-CSV-MA
load('log_BVAR_CSV_MA.mat', 'TC_L')
temp3(:,6)=median(TC_L,2); % this has weird values toward end of sample
% load in log_BVAR-CSV-t-MA
load('log_BVAR_CSV_t_MA.mat', 'TC_L')
temp3(:,7)=median(TC_L,2);
clear TC_L
end

C1=corr(temp);
C2=corr(temp2);
C3=corr(temp3);

X=table(dates,temp(:,1),temp(:,2),temp(:,3),temp(:,4),temp(:,5),temp(:,6),temp(:,7));
if LOG==0
writetable(X,'TC_S.csv')
elseif LOG==1
writetable(X,'TC_S_log.csv')
end

X=table(dates,temp2(:,1),temp2(:,2),temp2(:,3),temp2(:,4),temp2(:,5),temp2(:,6),temp2(:,7));
if LOG==0
writetable(X,'TC_M.csv')
elseif LOG==1
writetable(X,'TC_M_log.csv')
end

X=table(dates,temp3(:,1),temp3(:,2),temp3(:,3),temp3(:,4),temp3(:,5),temp3(:,6),temp3(:,7));
if LOG==0
writetable(X,'TC_L.csv')
elseif LOG==1
writetable(X,'TC_L_log.csv')
end

mycolor=[0.5 0.01 0.9];

figure(2)
subplot(3,1,1)
plot(temp(:,1),'b','LineWidth',1.5)
hold on,
plot(temp(:,2),'r','LineWidth',1.5)
plot(temp(:,3),'m','LineWidth',1.5)
plot(temp(:,4),'k','LineWidth',1.5)
plot(temp(:,5),'g','LineWidth',1.5)
plot(temp(:,6),'c','LineWidth',1.5)
plot(temp(:,7),'color',mycolor,'LineWidth',1.5)
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
plot(temp2(:,6),'c','LineWidth',1.5)
plot(temp2(:,7),'color',mycolor,'LineWidth',1.5)
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
plot(temp3(:,6),'c','LineWidth',1.5)
plot(temp3(:,7),'color',mycolor,'LineWidth',1.5)
ylabel('Long-term')
axis tight
grid on
legend('BVAR-t','BVAR-CSV', 'BVAR-MA(1)',...
    'BVAR-t-CSV','BVAR-t-MA(1)',' BVAR-CSV-MA(1)','BVAR-t-CSV-MA(1)','Location','NorthWest')

%
if LOG==0
load BVAR_CSV_t_MA.mat
elseif LOG==1
load log_BVAR_CSV_t_MA.mat
end
qq=[0.16 0.5 0.84];
temp=quantile(TC_S,qq,2);
temp2=quantile(TC_M,qq,2);
temp3=quantile(TC_L,qq,2);

X=table(dates,temp(:,1),temp(:,2),temp(:,3));
if LOG==0
writetable(X,'Short.csv')
elseif LOG==1
writetable(X,'Short_log.csv')
end
X=table(dates,temp2(:,1),temp2(:,2),temp2(:,3));
if LOG==0
writetable(X,'Medium.csv')
elseif LOG==1
writetable(X,'Medium_log.csv')
end
X=table(dates,temp3(:,1),temp3(:,2),temp3(:,3));
if LOG==0
writetable(X,'Long.csv')
elseif LOG==1
writetable(X,'Long_log.csv')
end

figure(3)
subplot(3,1,1)
plot(temp(:,2),'r-','LineWidth',1.5)
hold on,
plot(temp(:,1),'r--')
plot(temp(:,3),'r--')
ylabel('Short-term')
axis tight
grid on
subplot(3,1,2)
plot(temp2(:,2),'k-','LineWidth',1.5)
hold on,
plot(temp2(:,1),'k--')
plot(temp2(:,3),'k--')
ylabel('Medium-term')
axis tight
grid on
subplot(3,1,3)
plot(temp3(:,2),'b-','LineWidth',1.5)
hold on,
plot(temp3(:,1),'b--')
plot(temp3(:,3),'b--')
ylabel('Long-term')
axis tight
grid on

nd_s=quantile(ND_S,qq,3);
nd_m=quantile(ND_M,qq,3);
nd_l=quantile(ND_L,qq,3);

nds_vxefa=squeeze(nd_s(:,1,:));
ndm_vxefa=squeeze(nd_m(:,1,:));
ndl_vxefa=squeeze(nd_l(:,1,:));
nds_vix=squeeze(nd_s(:,2,:));
ndm_vix=squeeze(nd_m(:,2,:));
ndl_vix=squeeze(nd_l(:,2,:));
nds_ovx=squeeze(nd_s(:,3,:));
ndm_ovx=squeeze(nd_m(:,3,:));
ndl_ovx=squeeze(nd_l(:,3,:));
nds_vxxle=squeeze(nd_s(:,4,:));
ndm_vxxle=squeeze(nd_m(:,4,:));
ndl_vxxle=squeeze(nd_l(:,4,:));
nds_gvx=squeeze(nd_s(:,5,:));
ndm_gvx=squeeze(nd_m(:,5,:));
ndl_gvx=squeeze(nd_l(:,5,:));
nds_vxslv=squeeze(nd_s(:,6,:));
ndm_vxslv=squeeze(nd_m(:,6,:));
ndl_vxslv=squeeze(nd_l(:,6,:));

X=table(dates,nds_vxefa(:,1),nds_vxefa(:,2),nds_vxefa(:,3));
if LOG==0
writetable(X,'vxefa_nds.csv')
elseif LOG==1
writetable(X,'vxefa_nds_log.csv')
end
X=table(dates,ndm_vxefa(:,1),ndm_vxefa(:,2),ndm_vxefa(:,3));
if LOG==0
writetable(X,'vxefa_ndm.csv')
elseif LOG==1
writetable(X,'vxefa_ndm_log.csv')
end
X=table(dates,ndl_vxefa(:,1),ndl_vxefa(:,2),ndl_vxefa(:,3));
if LOG==0
writetable(X,'vxefa_ndl.csv')
elseif LOG==1
writetable(X,'vxefa_ndl_log.csv')    
end

X=table(dates,nds_vix(:,1),nds_vix(:,2),nds_vix(:,3));
if LOG==0
writetable(X,'vix_nds.csv')
elseif LOG==1
writetable(X,'vix_nds_log.csv')   
end
X=table(dates,ndm_vix(:,1),ndm_vix(:,2),ndm_vix(:,3));
if LOG==0
writetable(X,'vix_ndm.csv')
elseif LOG==1
writetable(X,'vix_ndm_log.csv')   
end
X=table(dates,ndl_vix(:,1),ndl_vix(:,2),ndl_vix(:,3));
if LOG==0
writetable(X,'vix_ndl.csv')
elseif LOG==1
writetable(X,'vix_ndl_log.csv')   
end

X=table(dates,nds_ovx(:,1),nds_ovx(:,2),nds_ovx(:,3));
if LOG==0
writetable(X,'ovx_nds.csv')
elseif LOG==1
writetable(X,'ovx_nds_log.csv')   
end
X=table(dates,ndm_ovx(:,1),ndm_ovx(:,2),ndm_ovx(:,3));
if LOG==0
writetable(X,'ovx_ndm.csv')
elseif LOG==1
writetable(X,'ovx_ndm_log.csv')   
end
X=table(dates,ndl_ovx(:,1),ndl_ovx(:,2),ndl_ovx(:,3));
if LOG==0
writetable(X,'ovx_ndl.csv')
elseif LOG==1
writetable(X,'ovx_ndl_log.csv')   
end

X=table(dates,nds_vxxle(:,1),nds_vxxle(:,2),nds_vxxle(:,3));
if LOG==0
writetable(X,'vxxle_nds.csv')
elseif LOG==1
writetable(X,'vxxle_nds_log.csv')    
end
X=table(dates,ndm_vxxle(:,1),ndm_vxxle(:,2),ndm_vxxle(:,3));
if LOG==0
writetable(X,'vxxle_ndm.csv')
elseif LOG==1
writetable(X,'vxxle_ndm_log.csv')    
end
X=table(dates,ndl_vxxle(:,1),ndl_vxxle(:,2),ndl_vxxle(:,3));
if LOG==0
writetable(X,'vxxle_ndl.csv')
elseif LOG==1
writetable(X,'vxxle_ndl_log.csv')    
end

X=table(dates,nds_gvx(:,1),nds_gvx(:,2),nds_gvx(:,3));
if LOG==0
writetable(X,'gvx_nds.csv')
elseif LOG==1
writetable(X,'gvx_nds_log.csv')
end
X=table(dates,ndm_gvx(:,1),ndm_gvx(:,2),ndm_gvx(:,3));
if LOG==0
writetable(X,'gvx_ndm.csv')
elseif LOG==1
writetable(X,'gvx_ndm_log.csv')
end
X=table(dates,ndl_gvx(:,1),ndl_gvx(:,2),ndl_gvx(:,3));
if LOG==0
writetable(X,'gvx_ndl.csv')
elseif LOG==1
writetable(X,'gvx_ndl_log.csv')
end

X=table(dates,nds_vxslv(:,1),nds_vxslv(:,2),nds_vxslv(:,3));
if LOG==0
writetable(X,'vxslv_nds.csv')
elseif LOG==1
writetable(X,'vxslv_nds_log.csv')   
end
X=table(dates,ndm_vxslv(:,1),ndm_vxslv(:,2),ndm_vxslv(:,3));
if LOG==0
writetable(X,'vxslv_ndm.csv')
elseif LOG==1
writetable(X,'vxslv_ndm_log.csv')   
end
X=table(dates,ndl_vxslv(:,1),ndl_vxslv(:,2),ndl_vxslv(:,3));
if LOG==0
writetable(X,'vxslv_ndl.csv')
elseif LOG==1
writetable(X,'vxslv_ndl_log.csv')   
end

figure(4)
subplot(3,2,1)
plot(nds_vxefa(:,2),'r-','Linewidth',1.5)
hold on,
plot(nds_vxefa(:,1),'r--')
plot(nds_vxefa(:,3),'r--')
axis tight
grid on
title('VXEFA')
subplot(3,2,2)
plot(nds_vix(:,2),'r-','Linewidth',1.5)
hold on,
plot(nds_vix(:,1),'r--')
plot(nds_vix(:,3),'r--')
axis tight
grid on
title('VIX')
subplot(3,2,3)
plot(nds_ovx(:,2),'r-','Linewidth',1.5)
hold on,
plot(nds_ovx(:,1),'r--')
plot(nds_ovx(:,3),'r--')
axis tight
grid on
title('OVX')
subplot(3,2,4)
plot(nds_vxxle(:,2),'r-','Linewidth',1.5)
hold on,
plot(nds_vxxle(:,1),'r--')
plot(nds_vxxle(:,3),'r--')
axis tight
grid on
title('VXXLE')
subplot(3,2,5)
plot(nds_gvx(:,2),'r-','Linewidth',1.5)
hold on,
plot(nds_gvx(:,1),'r--')
plot(nds_gvx(:,3),'r--')
axis tight
grid on
title('GVX')
subplot(3,2,6)
plot(nds_vxslv(:,2),'r-','Linewidth',1.5)
hold on,
plot(nds_vxslv(:,1),'r--')
plot(nds_vxslv(:,3),'r--')
axis tight
grid on
title('VXSLV')

figure(5)
subplot(3,2,1)
plot(ndm_vxefa(:,2),'k-','Linewidth',1.5)
hold on,
plot(ndm_vxefa(:,1),'k--')
plot(ndm_vxefa(:,3),'k--')
axis tight
grid on
title('VXEFA')
subplot(3,2,2)
plot(ndm_vix(:,2),'k-','Linewidth',1.5)
hold on,
plot(ndm_vix(:,1),'k--')
plot(ndm_vix(:,3),'k--')
axis tight
grid on
title('VIX')
subplot(3,2,3)
plot(ndm_ovx(:,2),'k-','Linewidth',1.5)
hold on,
plot(ndm_ovx(:,1),'k--')
plot(ndm_ovx(:,3),'k--')
axis tight
grid on
title('OVX')
subplot(3,2,4)
plot(ndm_vxxle(:,2),'k-','Linewidth',1.5)
hold on,
plot(ndm_vxxle(:,1),'k--')
plot(ndm_vxxle(:,3),'k--')
axis tight
grid on
title('VXXLE')
subplot(3,2,5)
plot(ndm_gvx(:,2),'k-','Linewidth',1.5)
hold on,
plot(ndm_gvx(:,1),'k--')
plot(ndm_gvx(:,3),'k--')
axis tight
grid on
title('GVX')
subplot(3,2,6)
plot(ndm_vxslv(:,2),'k-','Linewidth',1.5)
hold on,
plot(ndm_vxslv(:,1),'k--')
plot(ndm_vxslv(:,3),'k--')
axis tight
grid on
title('VXSLV')

figure(6)
subplot(3,2,1)
plot(ndl_vxefa(:,2),'b-','Linewidth',1.5)
hold on,
plot(ndl_vxefa(:,1),'b--')
plot(ndl_vxefa(:,3),'b--')
axis tight
grid on
title('VXEFA')
subplot(3,2,2)
plot(ndl_vix(:,2),'b-','Linewidth',1.5)
hold on,
plot(ndl_vix(:,1),'b--')
plot(ndl_vix(:,3),'b--')
axis tight
grid on
title('VIX')
subplot(3,2,3)
plot(ndl_ovx(:,2),'b-','Linewidth',1.5)
hold on,
plot(ndl_ovx(:,1),'b--')
plot(ndl_ovx(:,3),'b--')
axis tight
grid on
title('OVX')
subplot(3,2,4)
plot(ndl_vxxle(:,2),'b-','Linewidth',1.5)
hold on,
plot(ndl_vxxle(:,1),'b--')
plot(ndl_vxxle(:,3),'b--')
axis tight
grid on
title('VXXLE')
subplot(3,2,5)
plot(ndl_gvx(:,2),'b-','Linewidth',1.5)
hold on,
plot(ndl_gvx(:,1),'b--')
plot(ndl_gvx(:,3),'b--')
axis tight
grid on
title('GVX')
subplot(3,2,6)
plot(ndl_vxslv(:,2),'b-','Linewidth',1.5)
hold on,
plot(ndl_vxslv(:,1),'b--')
plot(ndl_vxslv(:,3),'b--')
axis tight
grid on
title('VXSLV')








