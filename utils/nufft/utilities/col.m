 function x = col(x)
%function x = col(x)
%	"colon" function
% x = x(:);
x = reshape(x,numel(x),1);