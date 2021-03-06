function [f] = get_normaldist(x,mu,sigma)

f = 1 ./ (sigma .* sqrt(2.*pi)) .* exp(-(x - mu).^2 ./ (2.*sigma.^2));
