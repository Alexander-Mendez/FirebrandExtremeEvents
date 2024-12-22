clear all

%% Inputs

%Number of integrations
k = 10;

%Number of realizations per integration
n = 20000;

%Initialize cell of landing locations
diastore = cell(1,k);
land = cell(1,k);

%%Firebrand parameters
rho_p = 513; %Firebrand density (kg/m^3)
mn = 5/1000; %Mean of diameters (m)
vr = 1/1000; %Standard deviation of diameters (m)
cd = 0.45; %Drag coefficient

%Log normal mean and variance
mu = log(mn^2/sqrt(vr+mn^2));
sigma = sqrt(log(vr/mn^2 + 1));

%%Flow parameters
rho_f = 1.204; %Ambient air density (kg/m^3)
U = 5; %Asymptotic velocity (m/s)
L = 25; %Boundary layer thickness (m)

%%Other parameters
g = 9.8; %Gravitational acceleration (m/s^2)
tspan = [0 42]; %Time span (s)
x0 = 0; %Initial horizontal position (m)
z0 = 50; %Initial vertical position (m)
v10 = 0; %Initial horizontal velocity (m/s)
v20 = 0; %Initial vertical velocity (m/s)

%Outputs vector of landing locations and the trajectories of the n
%firebrands from the k integration.

%% Integration loop
for i = 1:k
    % Log normally distributed diameters
    dia = lognrnd(mu,sigma,1,n);
    %Remove any negative diameters
    dia(dia<0)=NaN;
    diastore{1,i} = dia;
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
    [t,dat] = ode45(@(t,y) NDdvdt(t,y,p,q,n,@tanhflow), tspanhat, dat0);
    
    %% Determine firebrand landing position
    %Use a linear approximation, taking location of firebrand before and
    %after it crosses the x-axis then creating a line from these two points.
    %The x-intercept of the line is the approximate landing location.
    
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

    landk = x1 -z1'.*(x2-x1)./(z2-z1)' ;
    land{1,i} = landk;
    
    %Give trajectories in dimensional units
    dat(:,1:4:end) = dat(:,1:4:end)*L;
    dat(:,2:4:end) = dat(:,2:4:end)*L;
    dat(:,3:4:end) = dat(:,3:4:end)*U;
    dat(:,4:4:end) = dat(:,4:4:end)*U;
end
dia =cell2mat(diastore);
land = cell2mat(land);
save('dia.mat','dia')
save('land.mat','land')
save('traj.mat','dat')