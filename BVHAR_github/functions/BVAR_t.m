% This script estimates the BVAR-t model
%
% See:
% Chan, J.C.C. (2020). Large Bayesian VARs: A flexible Kronecker error 
% covariance structure, Journal of Business and Economic Statistics, 
% 38(1), 68-79.

%% prior
S0 = eye(n); nu0 = n+3;
construct_prior_A_network;
nuub = 100; %% upperbound for nu

%% construct X
X = zeros(T,n*p); 
for i=1:1
    X(:,(i-1)*n+1:i*n) = tmpY(p-i+1:end-i,:);
end
for i=5:5
    X(:,(i-1)*n+1:i*n) = tmpY(p-i+1:end-i,:);
end
for i=p:p
    X(:,(i-1)*n+1:i*n) = tmpY(p-i+1:end-i,:);
end
ind = find(sum(X,1)==0) ;
 X(:,ind) = [] ;
X = [ones(T,1) X];
   
%% initialize for storage
store_Sig = zeros(n,n,nsims); 
store_A = zeros(k,n,nsims);
store_nu = zeros(nsims,1); 
store_lam = zeros(nsims,T);
nugrid = linspace(2,nuub,700)';
store_pnu = zeros(700,1);

%% initialize the chain
nu = 5;
lam = 1./gamrnd(nu/2,2/nu,T,1);
countnu = 0;

%% MCMC starts here
randn('seed',sum(clock*100)); rand('seed',sum(clock*1000));
%disp('Starting MCMC for BVAR-t.... ');
start_time = clock;

for isim = 1:nsims + burnin
  
        %% sample Sig and A    
    iOm = sparse(1:T,1:T,1./lam);
    XiOm = X'*iOm;
    KA = sparse(1:k,1:k,1./VA0) + XiOm*X;
    Ahat = KA\(sparse(1:k,1:k,VA0)\A0 + XiOm*shortY);
    Shat = S0 + A0'*sparse(1:k,1:k,1./VA0)*A0 + shortY'*iOm*shortY ...
        - Ahat'*KA*Ahat;
    Shat = (Shat+Shat')/2; % adjust for rounding errors
    Sig = iwishrnd(Shat,nu0+T);
    CSig = chol(Sig,'lower');
    A = Ahat + (chol(KA,'lower')'\randn(k,n))*CSig'; 
    
        %% sample lam
    U = shortY - X*A;
    tmp = (U/CSig');
    s2 = sum(tmp.^2,2);
    lam = 1./gamrnd((n+nu)/2,2./(s2+nu));
    
        %% sample nu
    [nu,flag,fnu] = sample_nu(lam,nu,nuub);
    countnu = countnu + flag;

    if isim > burnin
        isave = isim - burnin; 
        store_A(:,:,isave) = A;
        store_Sig(:,:,isave) = Sig;        
        store_nu(isave,:) = nu;
        store_lam(isave,:) = lam';
        
        % scale sigma here as input to get_GIRF as per the paper. Note we
        % scale Sig by lam_{T} as Omega is a TxT diagonal matrix whose
        % elements are lam_{t} for t \in {1, 2, ..., T} This accounts for
        % the t-distribution in the theoretical quantities but also
        % correctly scales everything.
        
        scale=lam(end,1);
        Sig2=scale*Sig;
       
        
        HO=100+1;
        [~,wold]=get_GIRF(A',Sig2,1,L,HO-1);
        [ttfc,tc1,tc2,tc3,~,~,~,~,~,~,~,~,~,~,~,ndc1,ndc2,ndc3,tnd1]=get_dynnet(wold,T,Sig2,0); 
        % efficient code to get frequency connectedness
        tcs(:,isave)=tc3; tcm(:,isave)=tc2; tcl(:,isave)=tc1;
        tfc(:,isave)=ttfc;
        ndcl(:,isave)=ndc1; ndcm(:,isave)=ndc2; ndcs(:,isave)=ndc3; tndc(:,isave)=tnd1;
        % compute the density of nu
        %tmpden = fnu(nugrid);
        %tmpden = exp(tmpden-max(tmpden));
        %tmpden = tmpden/(sum(tmpden)*(nugrid(2)-nugrid(1)));
        %store_pnu = store_pnu + tmpden; 
    end
    
    %if ( mod(isim, 5000) ==0 )
    %    disp(  [ num2str(isim) ' loops... ' ] )
    %end 
    
end

%disp( ['MCMC takes '  num2str( etime( clock, start_time) ) ' seconds' ] );
%disp(' ' );

%A_mean = store_A/nsims;
%Sig_mean = store_Sig/nsims;
%nu_mean = mean(store_nu)';
%pnu_mean = store_pnu/nsims;

%figure;
%colormap('hsv');
%imagesc(A_mean);
%colorbar;
%box off;
%title('Heat map of the VAR coefficients');    

%figure; 
%plot(nugrid,pnu_mean); 
%title('Posterior density of nu');
%box off; xlim([0 30]);

if cp_ml
    ml_BVAR_t;
end