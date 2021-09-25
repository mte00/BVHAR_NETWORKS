function [timecon,fev]=get_timeconnect(N,HO,irf)

fev=vardecomp(N,HO,irf);
fev=fev(:,:,end);
FF=sum(fev);
FF=sum(FF);
timecon=100*(1-trace(fev)/FF);

end