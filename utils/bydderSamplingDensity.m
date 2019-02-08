function d = bydderSamplingDensity(nufftStructure,w,iter)
% Sampling density compensation as described and implemented in:
% Bydder et al, Evaluation of optimal density weighting for regridding,
% 2007
% https://www.sciencedirect.com/science/article/pii/S0730725X06002943

% Gather sampling locations from nufftStructure (only works if not created
% with 'table' argument, i.e. it has all sampling locations precalculated)
H = abs(nufftStructure.p);

if (nargin < 2 || isempty(w))
    w = 0.5 * max(abs(H(:)));
end
if (nargin <3 || isempty(iter))
    iter = 10;
end

[m, n]=size(H);
d=1./(H*(H'*ones(m,1)));
g=1./d+w^2;
r=H*(ones(n,1)-H'*d);
s=r./g;
new=r'*s;

for i=1:iter
    q=H*(H'*s)+w^2*s;
    alpha=new/(s'*q);
    d=d+alpha*s;
    r=r-alpha*q;
    q=r./g;
    old=new;
    new=r'*q;
    s=q+(new/old)*s;

end

end
