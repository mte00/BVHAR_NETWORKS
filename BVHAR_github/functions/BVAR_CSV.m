% This script estimates the BVAR-CSV model
%
% See:
% Chan, J.C.C. (2020). Large Bayesian VARs: A flexible Kronecker error 
% covariance structure, Journal of Business and Economic Statistics, 
% 38(1), 68-79.

%% prior
S0 = eye(n); nu0 = n+3;
construct_prior_A_network;
nuh0 = 5; Sh0 = .01*(nuh0-1); 
rho0 = .9; Vrho = .2^2;

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
store_Sig = zeros(n,n); 
store_A = zeros(k,n);
store_h = zeros(nsims,T);
store_theta = zeros(nsims,2);

%% initialize the chain
h = zeros(T,1);
rho = .8;
sigh2 = .1;
Hrho = speye(T) - rho*sparse(2:T,1:(T-1),ones(1,T-1),T,T);
counth = 0; countrho = 0;

%% MCMC starts here
randn('seed',sum(clock*100)); rand('seed',sum(clock*1000));
%disp('Starting MCMC for BVAR-CSV.... ');
start_time = clock;

for isim = 1:nsims + burnin
  
        %% sample Sig and A    
    iOh = sparse(1:T,1:T,exp(-h));
    XiOh = X'*iOh;
    KA = sparse(1:k,1:k,1./VA0) + XiOh*X;
    Ahat = KA\(sparse(1:k,1:k,VA0)\A0 + XiOh*shortY);
    Shat = S0 + A0'*sparse(1:k,1:k,1./VA0)*A0 + shortY'*iOh*shortY ...
        - Ahat'*KA*Ahat;
    Shat = (Shat+Shat')/2; % adjust for rounding errors
    Sig = iwishrnd(Shat,nu0+T);    
    CSig = chol(Sig,'lower');
    A = Ahat + (chol(KA,'lower')'\randn(k,n))*CSig'; 
    
        %% sample h
    U = shortY - X*A;
    tmp = (U/CSig');
    s2 = sum(tmp.^2,2);
    [h flag] = sample_h(s2,rho,sigh2,h,n);
    counth = counth + flag;
    
        %% sample sigh2
    eh = [h(1)*sqrt(1-rho^2);  h(2:end)-rho*h(1:end-1)];    
    sigh2 = 1/gamrnd(nuh0+T/2,1/(Sh0 + sum(eh.^2)/2));

        %% sample rho
    Krho = 1/Vrho + sum(h(1:T-1).^2)/sigh2;
    rhohat = Krho\(rho0/Vrho + h(1:T-1)'*h(2:T)/sigh2);
    rhoc = rhohat + sqrt(Krho)'\randn;
    grho = @(x) -.5*log(sigh2./(1-x.^2))-.5*(1-x.^2)/sigh2*h(1)^2;
    if abs(rhoc)<.9999
        alpMH = exp(grho(rhoc)-grho(rho));
        if alpMH>rand
            rho = rhoc;
            countrho = countrho+1;
            Hrho = speye(T) - rho*sparse(2:T,1:(T-1),ones(1,T-1),T,T);
        end
    end   
    
    if isim > burnin
        isave = isim - burnin; 
        store_A = store_A + A;
        store_Sig = store_Sig + Sig;
        store_h(isave,:) = h';
        store_theta(isave,:) = [rho sigh2];
        
        % scale sigma here as input to get_GIRF as per the paper. Note we
        % scale Sig by exp(h{T}) as Omega is a TxT diagonal matrix whose
        % elements are exp(h{t}) for t \in {1, 2, ..., T} This accounts for
        % the stochastic volatility in the theoretical quantities but also
        % correctly scales everything.
        
        scale=exp(h(end,1));
        Sig2=scale*Sig;
        
        HO=100+1;
        [~,wold]=get_GIRF(A',Sig2,1,L,HO-1);
        [ttfc,tc1,tc2,tc3,~,~,~,~,~,~,~,~,~,~,~,ndc1,ndc2,ndc3,tnd1]=get_dynnet(wold,T,Sig2,0); 
        % efficient code to get frequency connectedness
        tcs(:,isave)=tc3; tcm(:,isave)=tc2; tcl(:,isave)=tc1;
        tfc(:,isave)=ttfc;
        ndcl(:,isave)=ndc1; ndcm(:,isave)=ndc2; ndcs(:,isave)=ndc3; tndc(:,isave)=tnd1;
        
    end
    
    %if ( mod(isim, 5000) ==0 )
    %    disp(  [ num2str(isim) ' loops... ' ] )
    %end 
    
end

%disp( ['MCMC takes '  num2str( etime( clock, start_time) ) ' seconds' ] );
%disp(' ' );

A_mean = store_A/nsims;
Sig_mean = store_Sig/nsims;
h_mean = mean(store_h)';
theta_mean = mean(store_theta)';
CSV_std_mean = mean(exp(store_h/2))';

%figure;
%colormap('hsv');
%imagesc(A_mean);
%colorbar;
%box off;
%title('Heat map of the VAR coefficients');    

%figure; 
%plot(T_id, CSV_std_mean); box off;
%title('Posterior mean of exp(h_t/2)');
%box off; xlim([T_id(1)-1 T_id(end)+1]);

if cp_ml
    ml_BVAR_CSV;
end
