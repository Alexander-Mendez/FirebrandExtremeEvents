function temp = unclanding(z,lambda)
n =1;
dia = exp(z);

%%Firebrand parameters
rho_p = 513; %Firebrand density (kg/m^3)
mn = 1.5/1000; %Mean of diameters (m)
vr = (5/10000)^2; %Standard deviation of diameters (m)
cd = 0.45; %Drag coefficient

%Log normal mean and variance
mu = log(mn^2/sqrt(vr+mn^2));
sigma = sqrt(log(vr/mn^2 + 1));

%%Flow parameters
rho_f = 1.204; %Ambient air density (kg/m^3)
U = 0.7; %Asymptotic velocity (m/s)
L = 0.05; %Boundary layer thickness (m)
T = L/U;
%%Other parameters
g = 9.8; %Gravitational acceleration (m/s^2)
tspan = [0 650]; %Time span (s)
x0 = 0; %Initial horizontal position (m)
z0 = 50; %Initial vertical position (m)
v10 = 0; %Initial horizontal velocity (m/s)
v20 = 0; %Initial vertical velocity (m/s)

%Calculate dependent parameters
A = (pi/4)*(dia).^2; %Cross-sectional area
V = (pi/6)*(dia).^3; %Volume
m = V.*rho_p; %Mass of firebrand
T = L/U;
tspanhat = tspan./T;

p =[A;m];
q = [rho_f;cd;g;L;U];

dat0 = zeros(4,n);
dat0(1,:) = x0./L;
dat0(2,:) = (z0)./L;
dat0(3,:) = v10;
dat0(4,:) = v20;
opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
[t,dat] = ode45(@(t,y) NDdvdt2(t,y,p,q,n,@logflow), tspanhat, dat0,opts);


zdat = dat(:,2:4:end);

%Find index of first negative value
[z2ind, col] = find(cumsum(zdat<0)==1);

%Use the row index of the last positive vertical position
z1ind = z2ind -1;

z2 = zdat(sub2ind(size(zdat),z2ind,col));
z1 = zdat(sub2ind(size(zdat),z1ind,col));

%Use columns corresponding to horizontal position
col = (4*col - 2)-1;

%Return vector of landing distances
x2 = dat(sub2ind(size(dat),z2ind,col))'.*L(1);
x1 = dat(sub2ind(size(dat),z1ind,col))'.*L(1);

land = x1 -z1'.*(x2-x1)./(z2-z1)' ;

I = ((z - mu)^2)/(2*sigma^2);

temp = I - lambda*land;

end
