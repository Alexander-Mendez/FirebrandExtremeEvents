function [u] = logflow(t,xi1,xi2)

%% Calculate u
if xi2 > 1
    u1 = (1/0.4)*log(xi2);
    u2 = zeros(size(xi1));
    u = [u1;u2];
else
    u1 = zeros(size(xi1));
    u2 = zeros(size(xi1));
    u = [u1;u2];
end