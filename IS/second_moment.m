function [ val ] = second_moment( theta, mu, sigma, X)
val =  mean((theta(2)/sigma)*exp(-0.5*((X - mu)./(sigma)).^2 + 0.5*((X - theta(1))./(theta(2))).^2));
%kappa_theta = mu * theta + sigma^2 * theta^2 / 2;
%val = mean(exp(kappa_theta - theta*X));
end