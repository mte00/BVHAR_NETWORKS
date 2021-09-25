% This script constructs the prior for the VAR coefficients when getting
% network measures with HAR lag structure. 
%
% The method still follows See:
% Chan, J.C.C. (2020). Large Bayesian VARs: A flexible Kronecker error 
% covariance structure, Journal of Business and Economic Statistics, 
% 38(1), 68-79.

A0 = zeros(k,n);
VA0 = zeros(k,1);
sig2 = zeros(n,1);
tmpY = [Y0(end-p+1:end,:); shortY];
c1 = .2^2; c2 = 100; % hyperparameters for slopes and intercepts
for i=1:n
    Z = [ones(T,1) tmpY(22:end-1,i) tmpY(5:end-18,i) tmpY(1:end-22,i)];
    tmpb = (Z'*Z)\(Z'*tmpY(23:end,i));
    sig2(i) = mean((tmpY(23:end,i)-Z*tmpb).^2);
end
for i=1:k
    l = ceil((i-1)/n); 
    idx = mod(i-1,n); % variable index
    if idx==0
        idx = n;
    end        
    if i==1 % intercept
        VA0(1) = c2;     
    else
        VA0(i) = c1/(l^2*sig2(idx));        
    end   
end