function temp = NDdvdt(t,y,p,q,n,fun)
%%Parameters that depend on realizations of diameter.
Ac = p(1,:);
m = p(2,:);
%%Constant parameters
rho_f = q(1);
Cn = q(2);
g = q(3);
L = q(4);
U = q(5);
%% Nondimensionalization constants
gamma = (rho_f*Ac*Cn*L)./(2*m);
ghat = (L.*g)./(U.^2);

y = reshape(y,[],n);

[u] = fun(t,y(1,:),y(2,:));

v = [y(3,:);y(4,:)];

%Create gravity vector
gvec = [zeros(1,n);ghat*ones(1,n)];

%System of ODEs
dxi = v;
dv = (gamma.*vecnorm(u-v).*(u-v) - gvec);

temp = [dxi(1,:);dxi(2,:);dv(1,:);dv(2,:)];
temp = temp(:);
end