function [TFC,TC1,TC2,TC3,WC1,WC2,WC3]=get_connectedness4(wo,TT,sig)
%**************************************************************************
% Michael Ellington 28/11/2018
% Here this is following Defintion 2.3 of Barunik and Krehlik (2018).
% Noting that in the definition of (f(w))jk and the weighting function
% \Gammaj(w) the respective denominator and numerator cancel out. 
% Thereby simplifying the computation.
%
% Inputs: wo is wold decomposition of MA coefficients [N x N x HO]
%
%         TT is the number of time series observations used in estimating
%         TVP VAR model. This helps us create frequency window. I set to
%         10% of original time series sample since we have 2000+ daily
%         observations. This saves time!!!
%
%         sig is a draw of the estimated posterior distribution of the
%         covariance matrix of residuals stemming from the TVP VAR at time
%         t. [N x N]
%
% Outputs: TFC is total frequency connectedness
%          TC1 is long-term connectedness
%          TC2 is med-term connectedness
%          TC3 is short-term connectedness
%          WC1 is within long-term connectedness
%          WC2 is within med-term connectedness
%          WC3 is within short-term connectedness
%          
%          
%**************************************************************************
% Define frequency window;
Tw=floor(TT/10); % number of observations is 2611, this gives us 261 omegas.
omeg=linspace(0,pi,Tw)'; % create equally spaced line from 0 to pi in 261 intervals

% Define bands
omeg2=pi./omeg; 
bandmat=[omeg,omeg2];
d1=bandmat(bandmat(:,2)>20); % long term equals (20,260+] days
d2=bandmat(bandmat(:,2)<=20 & bandmat(:,2)>5); % medium term equals (5,20] days
d3=bandmat(bandmat(:,2)<=5); % short term equals [1,5] days
% d1 and d2 are not strictly open, since double counting when summing over
% frequencies will violate condition in line 155.

d1=length(d1); d2=length(d2); d3=length(d3);

% If omeg is frequency window, we decompose into three components:
% To get long-term we sum from 1:length(d1)
% To get med-term we sum from length(d1)+1:length(omeg)-length(d2)
% To get short-term we sum from
% length(d1)+length(d2)+1:end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here instead of frequency band, we are looking at (w={0,...,pi}) Which,
% due to symmetry of integral (-pi, pi) is the whole window...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HO=length(wo); % 3rd dimension of wold coefficients (i.e. 101)
N=size(wo,1);  % Number of variables (and shocks) in model
i=sqrt(-1);

% First get Omega = \sum_{w} \Psi(w)Sigma\Psi'(w) 
Omeg=zeros(N,N,length(omeg));
for w = 1:length(omeg)
    GI=zeros(N,N);
    for nn=1:HO
       GI=GI+wo(:,:,nn)*exp(-i*nn*omeg(w)); 
       %GI=GI+wo*exp(-i*nn*omeg(w)); 
    end
    Omeg(:,:,w)=GI*sig*GI'; 
end
Omeg=sum(real(Omeg),3); % This is overall omega (spec density over whole frequency window)
% The diagonal elements of this matrix are (Omeg)j,j that we use in line
% 96.


% First define numerator over all frequency windows
FC=zeros(N,N,length(omeg)); %FC1=FC;
for w=1:length(omeg)
    GI=zeros(N,N); % defining this as NaN produces equivalent value
   for nn=1:HO
        GI=GI+wo(:,:,nn)*exp(-i*nn*omeg(w)); % Takes fourier transform of wold at each h in 0,...,100 and sums them.
        %GI=GI+wo*exp(-i*nn*omeg(w)); % Takes fourier transform of wold at each h in 0,...,100 and sums them.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Here exp(-i*h*w) is slightly different to definition in Barunik
        % and Krehlik (2018) where they use exp(-{2*pi*i*w}/H), with
        % w={aH/2*pi, ... , bH/2*pi) with a<b defined as interval on real
        % line and each w rounded to lowest integer.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   end
   GI=GI*sig; % product of GI(w) and Sig at given frequency
   PS=zeros(N,N); PP=zeros(N,N);
   for k=1:N % following Defintion 2.3 in Barunik and Krehlik (2018)
       for j=1:N % jth row and kth column
        PS(j,k)=(abs(GI(j,k)))^2;
        PP(j,k)=PS(j,k)/(Omeg(j,j)*sig(k,k)); % Use Omeg as defined earlier and take ratio.
       end
   end
   
FC(:,:,w)=PP;
end

PP1=sum(FC,3);
for w=1:length(omeg)
for j=1:N
   FC(j,:,w)=FC(j,:,w)./sum(PP1(j,:)); % normalise theta by row sum of theta_inf 
end
end

thetainf=sum(FC,3);

% Now get long-term
temp1=sum(FC(:,:,1:d1),3); % theta_{d1} summed over (0 to 0.1576)
for j=1:N
    temp1(j,:)=temp1(j,:)./sum(thetainf(j,:));
end

%WC1=100*(1-trace(temp1)/sum(temp1,'all'));
%TC1=WC1*(sum(temp1,'all')/sum(thetainf,'all'));

WC1=100*(1-trace(temp1)/sum(sum(temp1)));
TC1=WC1*(sum(sum(temp1))/sum(sum(thetainf)));

% Now get med-term % theta_{d2} summed over (0.1576 to 0.6283)
temp1=sum(FC(:,:,d1+1:length(omeg)-d3),3);
for j=1:N
    temp1(j,:)=temp1(j,:)./sum(thetainf(j,:));
end

%WC2=100*(1-trace(temp1)/sum(temp1,'all'));
%TC2=WC2*(sum(temp1,'all')/sum(thetainf,'all'));
WC2=100*(1-trace(temp1)/sum(sum(temp1)));
TC2=WC2*(sum(sum(temp1))/sum(sum(thetainf)));

% Now get short-term % theta_{d3} summed over (0.6283 to \pi)
temp1=sum(FC(:,:,d1+d2+1:end),3);
for j=1:N
    temp1(j,:)=temp1(j,:)./sum(thetainf(j,:));
end

%WC3=100*(1-trace(temp1)/sum(temp1,'all'));
%TC3=WC3*(sum(temp1,'all')/sum(thetainf,'all'));
WC3=100*(1-trace(temp1)/sum(sum(temp1)));
TC3=WC3*(sum(sum(temp1))/sum(sum(thetainf)));
% Total Frequency Connect
tfc=TC1+TC2+TC3;

temp1=sum(FC,3);
   
   for j=1:N
       temp1(j,:)=temp1(j,:)./sum(temp1(j,:));
   end
% Total Frequency Connect: should be same as line 142.   
%tfc1=100*(sum(temp1,'all')/sum(thetainf,'all')-trace(temp1)/sum(thetainf,'all'));
tfc1=100*(sum(sum(temp1))/sum(sum(thetainf))-trace(temp1)/sum(sum(thetainf)));

if abs(tfc-tfc1)<eps+1.0e-10 % If these are equal to approximately 10d.p then allow function to go through.
   %TFC=tfc1;
   TFC=tfc;
else
    error('DIFFERENCE EXCEEDS ALLOWANCE OF eps+1.0e-10')
end

end