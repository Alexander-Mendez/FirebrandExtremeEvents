clear all
load('dia.mat')
load('land.mat')

%Define parameters of firebrand diameters
N_init = length(dia);
mu=5/1000;sigma=1/1000;
deltaL = 0.5;
iter = 1:16;
pdf_est = zeros(2,length(iter));
var_est = zeros(2,length(iter));
num_iter = zeros(2,length(iter));
pdf_MC = zeros(2,length(iter));
var_MC = zeros(2,length(iter));
theta_vec = zeros(2,length(iter));
num_dia = zeros(2,length(iter));


start = 39;
%pdf_est(1,:) = iter;
pdf_est(1,:) = iter+start;
var_est(1,:) = iter+start;
num_iter(1,:) = iter+start;
pdf_MC(1,:) = iter+start;
var_MC(1,:) = iter+start;
theta_vec(1,:) = iter+start;
num_dia(1,:) = iter+start;

for j = iter
Lstar = j + start;
%Lstar = j;
ind = abs(land-Lstar)<deltaL;
dia_sample = dia(ind);
L_sample = land(ind);
options = optimset('MaxFunEvals',1000000000);
theta_vm = fminsearch(@(theta) second_moment(theta,mu,sigma,dia_sample),[mu,sigma],options);
theta_vec(2,j) = theta_vm;
mu_is = mu + sigma^2*theta_vm;
n = 10^1;

dia_is = mu_is + sigma*randn(1,n);
dia_is = dia_is(dia_is>0);
n = length(dia_is);
num_dia(2,j) = n;

%Get samples of landing distributions in interval

%% 
%%Firebrand parameters
rho_p = 513; %Firebrand density (kg/m^3)
cd = 0.45; %Drag coefficient

%%Flow parameters
rho_f = 1.204; %Ambient air density (kg/m^3)
U = 5; %Asymptotic velocity (m/s)
L = 25; %Boundary layer thickness (m)

%%Other parameters
g = 9.8; %Gravitational acceleration (m/s^2)
tspan = [0 896]; %Time span (s)
x0 = 0; %Initial horizontal position (m)
z0 = 50; %Initial vertical position (m)
v10 = 0; %Initial horizontal velocity (m/s)
v20 = 0; %Initial vertical velocity (m/s)

A = (pi/4)*(dia_is).^2; %Cross-sectional area
V = (pi/6)*(dia_is).^3; %Volume
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
[t,dat] = ode45(@(t,y) NDdvdt(t,y,p,q,n,@tanhflow), tspanhat, dat0);

%Consider only columns giving vertical position
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

land_is = x1 -z1'.*(x2-x1)./(z2-z1)' ;

    
%Give trajectories in dimensional units
dat(:,1:4:end) = dat(:,1:4:end)*L;
dat(:,2:4:end) = dat(:,2:4:end)*L;
dat(:,3:4:end) = dat(:,3:4:end)*U;
dat(:,4:4:end) = dat(:,4:4:end)*U;

traject = [dat(:,1:4:end) dat(:,2:4:end)];

%%
results_IS = exp(mu*theta_vm + sigma^2*theta_vm^2/2 - theta_vm*dia_is).* (abs(land_is-Lstar)<deltaL);
results_MC = (abs(land-Lstar)<deltaL);

ell_hat_IS = mean(results_IS);
ell_hat_MC = mean(results_MC);

P_L = ell_hat_IS/(2*deltaL);
P_MC = ell_hat_MC/(2*deltaL);

pdf_est(2,j) = P_L;
pdf_MC(2,j) = P_MC;

var_est(2,j) = mean((results_IS - ell_hat_IS).^2);
var_MC(2,j) = mean((results_MC - ell_hat_MC).^2);
num_iter(2,j) = sum(abs(land_is-Lstar)<deltaL);
end

save('var_est.mat','var_est')
save('var_MC.mat','var_MC')
save('pdf_est.mat','pdf_est')
save('pdf_MC.mat','pdf_MC')
save('num_iter.mat','num_iter')
save('theta_vec.mat','theta_vec')
save('num_dia.mat','num_dia')
%save('raw.mat','raw')
save('traject.mat','traject')


%re_hat_IS = std(results_IS) / (sqrt(n)* ell_hat_IS);


