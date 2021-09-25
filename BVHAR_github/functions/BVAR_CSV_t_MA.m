% This script estimates the BVAR-CSV-t-MA model
%
% See:
% Chan, J.C.C. (2020). Large Bayesian VARs: A flexible Kronecker error 
% covariance structure, Journal of Business and Economic Statistics, 
% 38(1), 68-79.

%% prior
psi0 = 0; Vpsi = 1;
lpri_psi = @(x) -.5*(x-psi0)^2/Vpsi -1e10*(x<-.99 || x>.99);
S0 = eye(n); nu0 = n+3;
nuh0 = 5; Sh0 = .01*(nuh0-1); 
rho0 = .9; Vrho = .2^2;
nuub = 100; %% upperbound for nu
construct_prior_A_network;

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

%% initialize for storeage    
store_Sig = zeros(n,n); 
store_A = zeros(k,n);
store_h = zeros(nsims,T);
store_lam = zeros(nsims,T);
store_theta = zeros(nsims,4);

options = optimset('Display', 'off', 'LargeScale','off') ;
ngrid = 300;
psigrid = linspace(-.1,.4,ngrid)';
store_ppsi = zeros(ngrid,1);
nugrid = linspace(2,nuub,ngrid)';
store_pnu = zeros(ngrid,1);

%% initialize the chain
h = zeros(T,1);
nu = 5;
psi = -.1;
rho = .8;
sigh2 = .1;
Hrho = speye(T) - rho*sparse(2:T,1:(T-1),ones(1,T-1),T,T);
Hpsi = speye(T) + psi*sparse(2:T,1:(T-1),ones(1,T-1),T,T);
count_h = 0; count_rho = 0; count_psi = 0; count_nu = 0;
psihat = psi;
lam = 1./gamrnd(nu/2,2/nu,T,1);

%% MCMC starts here
randn('seed',sum(clock*100)); rand('seed',sum(clock*1000));
%disp('Starting MCMC for BVAR-CSV-t-MA.... ');
start_time = clock;

for isim = 1:nsims + burnin  
        % sample Sig and A   
    Xtld = Hpsi\X;    
    Ytld = Hpsi\shortY;
    iO_h_lam_psi = sparse(1:T,1:T,[1/(1+psi^2)*exp(-h(1));exp(-h(2:end))]./lam);
    XiO = Xtld'*iO_h_lam_psi;
    KA = sparse(1:k,1:k,1./VA0) + XiO*Xtld;
    Ahat = KA\(sparse(1:k,1:k,VA0)\A0 + XiO*Ytld);
    Shat = S0 + A0'*sparse(1:k,1:k,1./VA0)*A0 + Ytld'*iO_h_lam_psi*Ytld ...
        - Ahat'*KA*Ahat;
    Shat = (Shat+Shat')/2; % adjust for rounding errors
    Sig = iwishrnd(Shat,nu0+T);    
    CSig = chol(Sig,'lower');
    A = Ahat + (chol(KA,'lower')'\randn(k,n))*CSig'; 
    
        % sample lam
    U = shortY - X*A;
    Utld = Hpsi\U;
    tmp = (Utld/CSig');
    s2_lam = sum(tmp.^2,2)./exp(h);
    lam = 1./gamrnd((n+nu)/2,2./(s2_lam+nu));   
    
        % sample h    
    s2_h = sum(tmp.^2,2)./lam;
    s2_h(1) = s2_h(1)/(1+psi^2);    
    [h, flag] = sample_h(s2_h,rho,sigh2,h,n);
    count_h = count_h + flag;    
    
        % sample sigh2
    eh = [h(1)*sqrt(1-rho^2);  h(2:end)-rho*h(1:end-1)];    
    sigh2 = 1/gamrnd(nuh0+T/2,1/(Sh0 + sum(eh.^2)/2));

        % sample rho
    Krho = 1/Vrho + sum(h(1:T-1).^2)/sigh2;
    rhohat = Krho\(rho0/Vrho + h(1:T-1)'*h(2:T)/sigh2);
    rhoc = rhohat + sqrt(Krho)'\randn;
    grho = @(x) -.5*log(sigh2./(1-x.^2))-.5*(1-x.^2)/sigh2*h(1)^2;
    if abs(rhoc)<.99
        alpMH = exp(grho(rhoc)-grho(rho));
        if alpMH>rand
            rho = rhoc;
            count_rho = count_rho+1;
            Hrho = speye(T) - rho*sparse(2:T,1:(T-1),ones(1,T-1),T,T);
        end
    end 
    
        % sample psi
    U_psi = U./repmat(sqrt(lam),1,n);
    lp_psi = @(x) llike_CSV_MA(x,U_psi,Sig,h) + lpri_psi(x);    
    if (mod(isim,100)==0) || isim == 1 %% get the Hessian every 100 iterations
        [psihat,fval,exitflag,output,grad,hess] ...
            = fminunc(@(x)-lp_psi(x),psihat,options); 
        [tmpCpsi, flag] = chol(hess,'lower');
        if flag == 0
            Kpsic = hess;
        else 
            Kpsic = 1/.05^2;
        end
    else
        psihat = fminbnd(@(x)-lp_psi(x),-.99,.99);
    end    
    psic = psihat + 1/sqrt(Kpsic)*randn;
    if abs(psic)<.99
        alpMH =  lp_psi(psic) - lp_psi(psi) + ...
            -.5*(psi-psihat)^2*Kpsic + .5*(psic-psihat)^2*Kpsic;
    else
        alpMH = -inf;
    end
    if alpMH > log(rand)
        psi = psic;
        count_psi = count_psi + 1;
        Hpsi = speye(T) + psi*sparse(2:T,1:(T-1),ones(1,T-1),T,T);
    end     
    
        % sample nu
    [nu,flag,fnu] = sample_nu(lam,nu,nuub);
    count_nu = count_nu + flag;   
    
    if isim > burnin
        isave = isim - burnin; 
        store_A = store_A + A;
        store_Sig = store_Sig + Sig;
        store_h(isave,:) = h';
        store_lam(isave,:) = lam';
        store_theta(isave,:) = [psi rho sigh2 nu];
        
        % scale sigma here as input to get_GIRF as per the paper. Note we
        % scale Sig by (exp(h{T})+psi^2exp(h{T-1})) as Omega is a TxT diagonal matrix whose
        % elements follow svol MA(1) process.
        
        scale=(exp(h(end,1))*lam(end,1)+psi^2*exp(h(end-1,1))*lam(end-1,1));
        Sig2=scale*Sig;
        
        HO=100+1;
        [~,wold]=get_GIRF(A',Sig2,1,L,HO-1);
        [ttfc,tc1,tc2,tc3,~,~,~,~,~,~,~,~,~,~,~,ndc1,ndc2,ndc3,tnd1]=get_dynnet(wold,T,Sig2,0); 
        % efficient code to get frequency connectedness
        tcs(:,isave)=tc3; tcm(:,isave)=tc2; tcl(:,isave)=tc1;
        tfc(:,isave)=ttfc;
        ndcl(:,isave)=ndc1; ndcm(:,isave)=ndc2; ndcs(:,isave)=ndc3; tndc(:,isave)=tnd1;
        
        % Uncomment this block to compute the posterior density of psi
        % It would add a couple of minutes to the estimation        
%         % compute the density of psi
%         tmpden = zeros(ngrid,1);
%         for ii=1:ngrid
%             tmpden(ii) = lp_psi(psigrid(ii)); 
%         end
%         tmpden = exp(tmpden-max(tmpden));
%         tmpden = tmpden/(sum(tmpden)*(psigrid(2)-psigrid(1)));
%         store_ppsi = store_ppsi + tmpden;        


        % compute the density of nu
        tmpden = fnu(nugrid);
        tmpden = exp(tmpden-max(tmpden));
        tmpden = tmpden/(sum(tmpden)*(nugrid(2)-nugrid(1)));
        store_pnu = store_pnu + tmpden; 
    end
    
%    if ( mod(isim, 5000) ==0 )
%        disp(  [ num2str(isim) ' loops... ' ] )
%    end 
    
end

%disp( ['MCMC takes '  num2str( etime( clock, start_time) ) ' seconds' ] );
%disp(' ' );

A_mean = store_A/nsims;
Sig_mean = store_Sig/nsims;
CSV_std_mean = mean(exp(store_h/2))';
theta_mean = mean(store_theta)';
ppsi_mean = store_ppsi/nsims;
pnu_mean = store_pnu/nsims;

%figure;
%colormap('hsv');
%imagesc(A_mean);
%colorbar;
%box off;
%title('Heat map of the VAR coefficients');    

%T_id = linspace(1959,2013.75,T)';
%figure; 
%plot(T_id, CSV_std_mean); box off;
%title('Posterior mean of exp(h_t/2)');
%box off; xlim([T_id(1)-1 T_id(end)+1]);

%figure; 
%plot(nugrid,pnu_mean); 
%title('Posterior density of nu');
%box off; xlim([0 30]);

%figure; 
%hist(store_theta(:,1),50);
%title('Histogram of posterior draws of psi');
%box off;

% figure; 
% plot(psigrid,ppsi_mean); box off;
% title('Posterior density of psi');
% box off;

if cp_ml
    ml_BVAR_CSV_t_MA;
end




