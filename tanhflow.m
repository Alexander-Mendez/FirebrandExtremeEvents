function [u] = tanhflow(t,xi1,xi2)

%% Calculate u
u1 = (1 - (tanh(xi2)-1).^2);
u2 = zeros(size(xi1));
u = [u1;u2];
end