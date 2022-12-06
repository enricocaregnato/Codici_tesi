
function filt_coeffs=hyperfilter(hypermode,hyper_coeffs,degs,parms)

if nargin < 4 
    parms.lambda=10^(-1.5);
    parms.mu=ones(size(hyper_coeffs));
end

if isempty(parms)
    parms.lambda=10^(-1.5);
    parms.mu=ones(size(hyper_coeffs));
end

if length(parms.mu) == 1, parms.mu=parms.mu*ones(size(hyper_coeffs)); end

if size(degs,2) > size(degs,1), degs=degs'; end

switch hypermode
    case 'filtered'
        h=@(x) (x >= 0 & x <= 1/2) + ((sin(pi*x)).^2).*(x > 0.5 & x < 1);
        L=max(degs);
        filt_coeffs=h(degs/L).*hyper_coeffs;
    case 'lasso'
        lambda=parms.lambda; % regularization parameter
        mu=parms.mu; % penalties parameters
        S=@(a,k) max(0,a-k)+min(0,a+k);
        L=length(degs);
        for l=1:L
            al=hyper_coeffs(l);
            kl=lambda*mu(l);
            filt_coeffs(l,1)=S(al,kl);
        end
    otherwise
        filt_coeffs=hyper_coeffs;
end
