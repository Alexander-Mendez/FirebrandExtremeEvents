clear all
%lambda = [0.01:0.001:0.08];
%lambda = linspace(0.0001,0.6,1000);
%lambda = linspace(0.1,0.315,1000);
%lambda = linspace(0.055,0.2655,1000);
%lambda = linspace(0.040,0.3181,1500);
%lambda = [0.23];
%lambda = linspace(0.0001,0.2099,100);

%lambda = linspace(0.2654,0.27,1e5);
%lambda = [lambda(1), lambda(5000), lambda(10000), lambda(15000), lambda(20000), lambda(25000)];
%lambda2 = linspace(0.055,0.2655,10);
%lambda = [lambda2 lambda];

%lambda = linspace(0.055,0.27,4.7e6);
%lambda = linspace(0.040,0.26,300);
%lambda = linspace(0.10,0.29,300); %Goes to 270
%lambda = linspace(0.07,0.15,300); %Goes to 163
%lambda = lambda(1:296);
%lambda = linspace(0.07,0.11,30); %Goes to 240
lambda = linspace(0.07,0.11,3000); %Goes to 240
%lambda = linspace(lambda(19),0.27,4500);
%lambda = lambda(19);
est = zeros(3,length(lambda)+1);
est(1,1) = -7.366;
mn = 1.5/1000; %Mean of diameters (m)
vr = (5/10000)^2; %Standard deviation of diameters (m)

%Log normal mean and variance
mu = log(mn^2/sqrt(vr+mn^2));
sigma = sqrt(log(vr/mn^2 + 1));

options = optimset('TolFun',1e-8,'TolX',1e-8,'MaxFunEvals',100000000,'MaxIter',10000);

outiter = zeros(1,length(lambda));
outfunc = zeros(1,length(lambda));
tStart = tic;
for i = 1:length(lambda)
    z0 = est(1,i);
    f = @(z) unclanding(z,lambda(i));
    [zstar,val,exitflag,output] = fminunc(f,z0,options);
    I = ((zstar - mu)^2)/(2*sigma^2);
    P_Lstar = (2*pi)^(-1/2)*(exp(-I)/sqrt(2*I));    
    Lstar = (I-val)/lambda(i);

    est(1,i+1) = zstar;
    est(2,i+1) = Lstar;
    est(3,i+1) = P_Lstar;
    outiter(1,i) = output.iterations;
    outfunc(1,i) = output.funcCount;
end
tEnd = toc(tStart);

%figure
%semilogy(truncest(2,:),truncest(3,:),'-o','linewidth',2)
%xlabel('Threshold $L^*$ (m)','interpreter','latex','fontsize',20)
%ylabel('Probability','interpreter','latex','fontsize',20)

%figure
%semilogy(truncest(2,1:end-1),diff(truncest(3,:))./(-diff(truncest(2,:))));
%xlim([58 72])
%xlabel('Threshold $L^*$ (m)','interpreter','latex','fontsize',20)
%ylabel('Probability','interpreter','latex','fontsize',20)

%figure
%plot(lambda(1:57),est(1,2:58),'LineWidth',2)
%xlabel('$\lambda$','interpreter','latex','fontsize',22)
%ylabel('$z^*$','interpreter','latex','fontsize',22)

%figure
%plot(lambda(1:57),exp(est(1,2:58)),'LineWidth',2)
%xlabel('$\lambda$','interpreter','latex','fontsize',22)
%ylabel('Diameter (m)','interpreter','latex','fontsize',22)

save('est.mat','est')
save('outiter.mat','outiter')
save('outfunc.mat','outfunc')
save('tEnd.mat','tEnd')

save('est176to239.mat','est')
save('outiter176to239.mat','outiter')
save('outfunc176to239.mat','outfunc')
save('tEnd176to239.mat','tEnd')