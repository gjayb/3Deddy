%intuition plot of trajectory ensemble, in psi
%ranges in psi of trajectory ensemble
%% intuition figure

load('E0.125sigma0eps0.01r0.1Sb.mat')
psitrSteady=mypsi(R,E,xtr,ytr,ztr);
load('E0.125sigma0.042eps0r0.1S.mat')
psitrStoch=mypsi(R,E,xtr,ytr,ztr);
load('E0.125sigma0.042eps0.01r0.1Sb.mat')
psitrBoth=mypsi(R,E,xtr,ytr,ztr);

figure
subplot(3,1,1)
plot(ttr,psitrSteady,'Color',[0.5 0.5 0.5])
hold on
plot(ttr,mean(psitrSteady,2),'k','LineWidth',2)
plot(ttr,mean(psitrSteady,2)+std(psitrSteady,0,2),'k-.','LineWidth',2)
plot(ttr,mean(psitrSteady,2)-std(psitrSteady,0,2),'k-.','LineWidth',2)
ylim([-1e-4 1.5e-3])
xlim([0 3000])
title('E=0.125, Steady \epsilon=0.01','fontsize',20)
ylabel('\psi','fontsize',20)
set(gca,'fontsize',20)
subplot(3,1,2)
plot(ttr,psitrStoch,'Color',[0.5 0.5 0.5])
hold on
plot(ttr,mean(psitrStoch,2),'k','LineWidth',2)
plot(ttr,mean(psitrStoch,2)+std(psitrStoch,0,2),'k-.','LineWidth',2)
plot(ttr,mean(psitrStoch,2)-std(psitrStoch,0,2),'k-.','LineWidth',2)
ylim([-1e-4 1.5e-3])
xlim([0 3000])
title('Stochastic, \kappa=10^{-6}','fontsize',20)
set(gca,'fontsize',20)
ylabel('\psi','fontsize',20)
subplot(3,1,3)
plot(ttr,psitrBoth,'Color',[0.5 0.5 0.5])
hold on
plot(ttr,mean(psitrBoth,2),'k','LineWidth',2)
plot(ttr,mean(psitrBoth,2)+std(psitrBoth,0,2),'k-.','LineWidth',2)
plot(ttr,mean(psitrBoth,2)-std(psitrBoth,0,2),'k-.','LineWidth',2)
ylim([-1e-4 1.5e-3])
xlim([0 3000])
title('Both perturbations, \kappa=10^{-6} and \epsilon=0.01','fontsize',20)
set(gca,'fontsize',20)
xlabel('time','fontsize',20)
ylabel('\psi','fontsize',20)
%% E=0.125
% load('E0.125sigma0.042eps0.01r0.9Sphere.mat')
% psitrBoth=mypsi(R,E,xtr,ytr,ztr);
% load('E0.125sigma0eps0.01r0.9Sphere.mat')
% psitrSteady=mypsi(R,E,xtr,ytr,ztr);
% load('E0.125sigma0.042eps0r0.9Sphere.mat')
% psitrStoch=mypsi(R,E,xtr,ytr,ztr);

load('E0.125sigma0.013eps0.01r0.1S.mat')
psitrBoth=mypsi(R,E,xtr,ytr,ztr);
load('E0.125sigma0eps0.01r0.1S.mat')
psitrSteady=mypsi(R,E,xtr,ytr,ztr);
load('E0.125sigma0eps0.08r0.1S.mat')
psitrSteady2=mypsi(R,E,xtr,ytr,ztr);
load('E0.125sigma0.042eps0r0.1S.mat')
psitrStoch=mypsi(R,E,xtr,ytr,ztr);

figure%(1)
subplot(4,1,1)
plot(ttr,psitrSteady,'Color',[0.5 0.5 0.5])
hold on
plot(ttr,mean(psitrSteady,2),'k','LineWidth',2)
plot(ttr,mean(psitrSteady,2)+std(psitrSteady,0,2),'k-.','LineWidth',2)
plot(ttr,mean(psitrSteady,2)-std(psitrSteady,0,2),'k-.','LineWidth',2)
ylim([0 2e-3])
title('E=0.125, Steady \epsilon=0.01','fontsize',14)
set(gca,'fontsize',14)
subplot(4,1,2)
plot(ttr,psitrSteady2,'Color',[0.5 0.5 0.5])
hold on
plot(ttr,mean(psitrSteady2,2),'k','LineWidth',2)
plot(ttr,mean(psitrSteady2,2)+std(psitrSteady2,0,2),'k-.','LineWidth',2)
plot(ttr,mean(psitrSteady2,2)-std(psitrSteady2,0,2),'k-.','LineWidth',2)
ylim([0 2e-3])
title('Steady \epsilon=0.08','fontsize',14)
set(gca,'fontsize',14)
subplot(4,1,3)
plot(ttr,psitrBoth,'Color',[0.5 0.5 0.5])
hold on
plot(ttr,mean(psitrBoth,2),'k','LineWidth',2)
plot(ttr,mean(psitrBoth,2)+std(psitrBoth,0,2),'k-.','LineWidth',2)
plot(ttr,mean(psitrBoth,2)-std(psitrBoth,0,2),'k-.','LineWidth',2)
ylim([0 2e-3])
title('Steady, \epsilon=0.01, and Stochastic, \kappa=10^{-7}','fontsize',14)
set(gca,'fontsize',14)
subplot(4,1,4)
plot(ttr,psitrStoch,'Color',[0.5 0.5 0.5])
hold on
plot(ttr,mean(psitrStoch,2),'k','LineWidth',2)
plot(ttr,mean(psitrStoch,2)+std(psitrStoch,0,2),'k-.','LineWidth',2)
plot(ttr,mean(psitrStoch,2)-std(psitrStoch,0,2),'k-.','LineWidth',2)
ylim([0 2e-3])
title('Stochastic, \kappa=10^{-6}','fontsize',14)
set(gca,'fontsize',14)
xlabel('time','fontsize',14)


%load('E0.125sigma0.021eps0.005r0.9Sphere.mat') %doesnt exist (yet)
%psitrHalf=mypsi(R,E,xtr,ytr,ztr);
load('E0.125sigma0.013eps0r0.9Sphere.mat')
psitrStochSmall=mypsi(R,E,xtr,ytr,ztr);

figure%(2)
plot(ttr,range(psitrSteady,2),'LineWidth',2)
hold all
plot(ttr,range(psitrStoch,2),'LineWidth',2)
plot(ttr,range(psitrBoth,2),'LineWidth',2)
plot(ttr,range(psitrStochSmall,2),'LineWidth',2)
legend('steady, \epsilon=0.01','stochastic, \kappa=10^{-6}','both', '\kappa=10^{-7}')
xlabel('time','fontsize',14)
set(gca,'fontsize',14)
ylabel('range in \Psi','fontsize',14)
title('Range of Ensembles in \Psi, E=0.125, outer')

figure
plot(ttr,std(psitrSteady,0,2),'LineWidth',2)
hold all
plot(ttr,std(psitrStoch,0,2),'LineWidth',2)
plot(ttr,std(psitrBoth,0,2),'LineWidth',2)
%plot(ttr,range(psitrStoch0,2),'LineWidth',2)
legend('steady, \epsilon=0.01','stochastic, \kappa=10^{-6}','both', '\kappa=10^{-6}')
xlabel('time','fontsize',14)
set(gca,'fontsize',14)
ylabel('std(\Psi)','fontsize',14)
title('Standard deviation of Ensembles in \Psi, E=0.125, outer')
%% E=0.02
% load('E0.02sigma0.042eps0.01r0.9Sphere.mat')
% psitrBoth=mypsi(R,E,xtr,ytr,ztr);
% load('E0.02sigma0eps0.01r0.9Sphere.mat')
% psitrSteady=mypsi(R,E,xtr,ytr,ztr);
% load('E0.02sigma0.042eps0r0.9Sphere.mat')
% psitrStoch=mypsi(R,E,xtr,ytr,ztr);

load('E0.02sigma0eps0.08r0.1S.mat')
psitrBoth=mypsi(R,E,xtr,ytr,ztr);
load('E0.02sigma0eps0.01r0.1S.mat')
psitrSteady=mypsi(R,E,xtr,ytr,ztr);
load('E0.02sigma0.042eps0r0.1S.mat')
psitrStoch=mypsi(R,E,xtr,ytr,ztr);

figure%(3)
subplot(3,1,1)
plot(ttr,psitrSteady,'Color',[0.5 0.5 0.5])
hold on
plot(ttr,mean(psitrSteady,2),'k','LineWidth',2)
plot(ttr,mean(psitrSteady,2)+std(psitrSteady,0,2),'k-.','LineWidth',2)
plot(ttr,mean(psitrSteady,2)-std(psitrSteady,0,2),'k-.','LineWidth',2)
ylim([0 2e-3])
title('E=0.02, Steady \epsilon=0.01','fontsize',14)
set(gca,'fontsize',14)
subplot(3,1,2)
plot(ttr,psitrBoth,'Color',[0.5 0.5 0.5])
hold on
plot(ttr,mean(psitrBoth,2),'k','LineWidth',2)
plot(ttr,mean(psitrBoth,2)+std(psitrBoth,0,2),'k-.','LineWidth',2)
plot(ttr,mean(psitrBoth,2)-std(psitrBoth,0,2),'k-.','LineWidth',2)
ylim([0 2e-3])
title('Steady, \epsilon=0.08','fontsize',14)
set(gca,'fontsize',14)
subplot(3,1,3)
plot(ttr,psitrStoch,'Color',[0.5 0.5 0.5])
hold on
plot(ttr,mean(psitrStoch,2),'k','LineWidth',2)
plot(ttr,mean(psitrStoch,2)+std(psitrStoch,0,2),'k-.','LineWidth',2)
plot(ttr,mean(psitrStoch,2)-std(psitrStoch,0,2),'k-.','LineWidth',2)
ylim([0 2e-3])
xlabel('time','fontsize',14)
set(gca,'fontsize',14)
title('Stochastic, \kappa=10^{-6}','fontsize',14)


%load('E0.125sigma0.021eps0.005r0.9Sphere.mat') %doesnt exist (yet)
%psitrHalf=mypsi(R,E,xtr,ytr,ztr);
load('E0.02sigma0.013eps0r0.9Sphere.mat')
psitrStochSmall=mypsi(R,E,xtr,ytr,ztr);

figure%(4)
plot(ttr,range(psitrSteady,2),'LineWidth',2)
hold all
plot(ttr,range(psitrStoch,2),'LineWidth',2)
plot(ttr,range(psitrBoth,2),'LineWidth',2)
plot(ttr,range(psitrStochSmall,2),'LineWidth',2)
legend('steady, \epsilon=0.01','stochastic, \kappa=10^{-6}','both', '\kappa=10^{-7}')
xlabel('time','fontsize',14)
set(gca,'fontsize',14)
ylabel('range in \Psi','fontsize',14)
title('Range of Ensembles in \Psi, E=0.02, outer')

%% E=0.125 inner
load('E0.125sigma0.042eps0.01r0.5Sphere.mat')
psitrBoth=mypsi(R,E,xtr,ytr,ztr);
load('E0.125sigma0eps0.01r0.5Sphere.mat')
psitrSteady=mypsi(R,E,xtr,ytr,ztr);
load('E0.125sigma0.042eps0r0.5Sphere.mat')
psitrStoch=mypsi(R,E,xtr,ytr,ztr);
load('E0.125sigma0.013eps0r0.5Sphere.mat')
psitrStochSmall=mypsi(R,E,xtr,ytr,ztr);

figure%(5)
plot(ttr,range(psitrSteady,2),'LineWidth',2)
hold all
plot(ttr,range(psitrStoch,2),'LineWidth',2)
plot(ttr,range(psitrBoth,2),'LineWidth',2)
plot(ttr,range(psitrStochSmall,2),'LineWidth',2)
legend('steady, \epsilon=0.01','stochastic, \kappa=10^{-6}','both', '\kappa=10^{-7}')
xlabel('time','fontsize',14)
set(gca,'fontsize',14)
ylabel('range in \Psi','fontsize',14)
title('Range of Ensembles in \Psi, E=0.125, inner')
%% E=0.02 inner
load('E0.02sigma0.042eps0.01r0.5Sphere.mat')
psitrBoth=mypsi(R,E,xtr,ytr,ztr);
load('E0.02sigma0eps0.01r0.5Sphere.mat')
psitrSteady=mypsi(R,E,xtr,ytr,ztr);
load('E0.02sigma0.042eps0r0.5Sphere.mat')
psitrStoch=mypsi(R,E,xtr,ytr,ztr);
load('E0.02sigma0.013eps0r0.5Sphere.mat')
psitrStochSmall=mypsi(R,E,xtr,ytr,ztr);

figure%(6)
plot(ttr,range(psitrSteady,2),'LineWidth',2)
hold all
plot(ttr,range(psitrStoch,2),'LineWidth',2)
plot(ttr,range(psitrBoth,2),'LineWidth',2)
plot(ttr,range(psitrStochSmall,2),'LineWidth',2)
legend('steady, \epsilon=0.01','stochastic, \kappa=10^{-6}','both', '\kappa=10^{-7}')
xlabel('time','fontsize',14)
set(gca,'fontsize',14)
ylabel('range in \Psi','fontsize',14)
title('Range of Ensembles in \Psi, E=0.02, inner')

%% E=0.125 outer various z

load('E0.125sigma0eps0.01r0.9Sphere.mat')
% pert=4; x0=-0.5; trajectorySlice
psitrSteady=mypsi(R,E,xtr,ytr,ztr);
load('E0.125sigma0eps0.01y0.9z0.1Sphere.mat')
%trajectorySlice
psitrZ1=mypsi(R,E,xtr,ytr,ztr);
load('E0.125sigma0eps0.01y0.9z0.3Sphere.mat')
%trajectorySlice
psitrZ3=mypsi(R,E,xtr,ytr,ztr);
load('E0.125sigma0eps0.01y0.9z0.7Sphere.mat')
%trajectorySlice
psitrZ7=mypsi(R,E,xtr,ytr,ztr);
load('E0.125sigma0eps0.01y0.9z0.85Sphere.mat')
%trajectorySlice
psitrZ85=mypsi(R,E,xtr,ytr,ztr);
load('E0.125sigma0eps0.01y0.9z0.95Sphere.mat')
%trajectorySlice
psitrZ95=mypsi(R,E,xtr,ytr,ztr);

figure
plot(ttr,range(psitrSteady,2),'LineWidth',2)
hold all
plot(ttr,range(psitrZ1,2),'LineWidth',2)
plot(ttr,range(psitrZ3,2),'LineWidth',2)
plot(ttr,range(psitrZ7,2),'LineWidth',2)
plot(ttr,range(psitrZ85,2),'LineWidth',2)
plot(ttr,range(psitrZ95,2),'LineWidth',2)
xlabel('time','fontsize',14)
set(gca,'fontsize',14)
ylabel('range in \Psi','fontsize',14)
title('Range of Ensembles in \Psi, E=0.125, y=0.9')

%% E=0.02 outer various z

load('E0.02sigma0eps0.01r0.9Sphere.mat')
%trajectorySlice
psitrSteady=mypsi(R,E,xtr,ytr,ztr);
load('E0.02sigma0eps0.01y0.9z0.1Sphere.mat')
%trajectorySlice
psitrZ1=mypsi(R,E,xtr,ytr,ztr);
load('E0.02sigma0eps0.01y0.9z0.3Sphere.mat')
%trajectorySlice
psitrZ3=mypsi(R,E,xtr,ytr,ztr);
load('E0.02sigma0eps0.01y0.9z0.7Sphere.mat')
%trajectorySlice
psitrZ7=mypsi(R,E,xtr,ytr,ztr);
load('E0.02sigma0eps0.01y0.9z0.85Sphere.mat')
%trajectorySlice
psitrZ85=mypsi(R,E,xtr,ytr,ztr);
load('E0.02sigma0eps0.01y0.9z0.95Sphere.mat')
%trajectorySlice
psitrZ95=mypsi(R,E,xtr,ytr,ztr);

figure
plot(ttr,range(psitrSteady,2),'LineWidth',2)
hold all
plot(ttr,range(psitrZ1,2),'LineWidth',2)
plot(ttr,range(psitrZ3,2),'LineWidth',2)
plot(ttr,range(psitrZ7,2),'LineWidth',2)
plot(ttr,range(psitrZ85,2),'LineWidth',2)
plot(ttr,range(psitrZ95,2),'LineWidth',2)
legend('steady, z=0.5','z=0.1','z=0.3','z=0.7','z=0.85','z=0.95')
xlabel('time','fontsize',14)
set(gca,'fontsize',14)
ylabel('range in \Psi','fontsize',14)
title('Range of Ensembles in \Psi, E=0.02, y=0.9')

%% E=0.125 central line (y=0) various z
pert=4; x0=-0.5;
load('E0.125sigma0eps0.01y0z0.1Sphere.mat')
%trajectorySlice
psitrZ1=mypsi(R,E,xtr,ytr,ztr);
load('E0.125sigma0eps0.01y0z0.3Sphere.mat')
%trajectorySlice
psitrZ3=mypsi(R,E,xtr,ytr,ztr);
load('E0.125sigma0eps0.01y0z0.7Sphere.mat')
%trajectorySlice
psitrZ7=mypsi(R,E,xtr,ytr,ztr);
load('E0.125sigma0eps0.01y0z0.9Sphere.mat')
%trajectorySlice
psitrZ9=mypsi(R,E,xtr,ytr,ztr);


figure
plot(ttr,range(psitrZ1,2),'LineWidth',2)
hold all
plot(ttr,range(psitrZ3,2),'LineWidth',2)
plot(ttr,range(psitrZ7,2),'LineWidth',2)
plot(ttr,range(psitrZ9,2),'LineWidth',2)
legend('z=0.1','z=0.3','z=0.7','z=0.9')
xlabel('time','fontsize',14)
set(gca,'fontsize',14)
ylabel('range in \Psi','fontsize',14)
title('Range of Ensembles in \Psi, E=0.125, y=0')

%% E=0.02 central line (y=0) various z
pert=4; x0=-0.5;
load('E0.02sigma0eps0.01y0z0.1Sphere.mat')
%trajectorySlice
psitrZ1=mypsi(R,E,xtr,ytr,ztr);
load('E0.02sigma0eps0.01y0z0.3Sphere.mat')
%trajectorySlice
psitrZ3=mypsi(R,E,xtr,ytr,ztr);
load('E0.02sigma0eps0.01y0z0.7Sphere.mat')
%trajectorySlice
psitrZ7=mypsi(R,E,xtr,ytr,ztr);
load('E0.02sigma0eps0.01y0z0.9Sphere.mat')
%trajectorySlice
psitrZ9=mypsi(R,E,xtr,ytr,ztr);


figure
plot(ttr,range(psitrZ1,2),'LineWidth',2)
hold all
plot(ttr,range(psitrZ3,2),'LineWidth',2)
plot(ttr,range(psitrZ7,2),'LineWidth',2)
plot(ttr,range(psitrZ9,2),'LineWidth',2)
legend('z=0.1','z=0.3','z=0.7','z=0.9')
xlabel('time','fontsize',14)
set(gca,'fontsize',14)
ylabel('range in \Psi','fontsize',14)
title('Range of Ensembles in \Psi, E=0.02, y=0')


%% make movie, 3D
pert=4; x0=-0.5;
%load('E0.125sigma0eps0.01y0z0.1Sphere.mat')
%trajectorySlice
load('E0.125sigma0.042eps0y0.05z0.1SphereV3.mat')
psitrZ1=mypsi(R,E,xtr,ytr,ztr);

v = VideoWriter('psiSpreadTraj0125CenterlineStochV3.avi');
v.FrameRate=10;
open(v)

    f1=figure('renderer','zbuffer');
for i=6:2:5000
    clear testframe
    plot3(xtr(i-5:i,:),ytr(i-5:i,:),ztr(i-5:i,:))
    hold on
    plot3(xtr(i,:),ytr(i,:),ztr(i,:),'o')
    axis([-1 1 -1 1 0 1])
    grid on

title(strcat('Trajectories, \epsilon=0.01, E=0.125, t=',num2str(i),' range(\psi)=',num2str(range(psitrZ1(i,:)))))


testframe=getframe(gcf);
writeVideo(v,getframe(gcf))
writeVideo(v,testframe)

 clf   
end
close(v)
close(f1)

%% make movie, 2D
pert=4; x0=-0.5;
%load('E0.125sigma0eps0.01y0z0.1Sphere.mat')
%load('E0.02sigma0eps0.01y0z0.1Sphere.mat')
%load('E0.125sigma0.042eps0y0.05z0.1SphereV3.mat')
load('E0.125sigma0eps0y0.05z0.1Sphere.mat')
trajectorySlice
psitrZ1=mypsi(R,E,xtr,ytr,ztr);
[thetatr,rtr]=cart2pol(xtr,ytr);

r1=-1:0.01:1;
z1=0:0.01:1;
[rr1,zz1]=meshgrid(r1,z1);
psi1=mypsi(R,E,rr1,zeros(size(rr1)),zz1);

v = VideoWriter('psiSpreadTraj0125planeCenterlineNopert.avi');
v.FrameRate=10;
open(v)

    f1=figure('renderer','zbuffer');
for i=6:2:5000
    clear testframe
    contour(rr1,zz1,psi1)
    hold on
    plot(rtr(i-5:i,:).*repmat(sign(xtr(i,:)),[6 1]),ztr(i-5:i,:))
    plot(rtr(i,:).*sign(xtr(i,:)),ztr(i,:),'o')
    axis([-1 1 0 1])
    grid on

title(strcat('Trajectories, \epsilon=0, E=0.125, t=',num2str(i),' range(\psi)=',num2str(range(psitrZ1(i,:)))))


testframe=getframe(gcf);
writeVideo(v,getframe(gcf))
writeVideo(v,testframe)

 clf   
end
close(v)
close(f1)

%% E=0.125 outer, starting at x,y=0, z=0.1
pert=4; x0=-0.5;
%load('E0.125sigma0.042eps0.01y0z0.1SphereV3.mat')
%psitrBoth1=mypsi(R,E,xtr,ytr,ztr);
%[~,rB]=cart2pol(xtr,ytr);
%trajectorySlice

%load('E0.125sigma0.042eps0y0z0.1SphereV2.mat')
%psitrStoch=mypsi(R,E,xtr,ytr,ztr);
%trajectorySlice
%load('E0.125sigma0.042eps0y0z0.1SphereV3.mat')
%psitrStoch1=mypsi(R,E,xtr,ytr,ztr);
%[~,rR]=cart2pol(xtr,ytr);
%trajectorySlice

%load('E0.125sigma0eps0.01y0z0.1Sphere.mat')
%psitrSteady1=mypsi(R,E,xtr,ytr,ztr);
%[~,rS]=cart2pol(xtr,ytr);
%trajectorySlice

load('E0.125sigma0.042eps0.01y0z0.1Sphere.mat')
psitrBoth1=mypsi(R,E,xtr,ytr,ztr);
load('E0.125sigma0.042eps0y0z0.1Sphere.mat')
psitrStoch1=mypsi(R,E,xtr,ytr,ztr);
load('E0.125sigma0eps0.01y0z0.1Sphere.mat')
psitrSteady1=mypsi(R,E,xtr,ytr,ztr);

figure
plot(ttr,range(psitrSteady1,2),'LineWidth',2)
hold all
plot(ttr,range(psitrStoch1,2),'LineWidth',2)
plot(ttr,range(psitrBoth1,2),'LineWidth',2)
%plot(ttr,range(psitrStoch0,2),'LineWidth',2)
legend('steady, \epsilon=0.01','stochastic, \kappa=10^{-6}','both', '\kappa=10^{-6}')
xlabel('time','fontsize',14)
set(gca,'fontsize',14)
ylabel('range in \Psi','fontsize',14)
title('Range of Ensembles in \Psi, E=0.125, outer')

figure
plot(ttr,std(psitrSteady1,0,2),'LineWidth',2)
hold all
plot(ttr,std(psitrStoch1,0,2),'LineWidth',2)
plot(ttr,std(psitrBoth1,0,2),'LineWidth',2)
%plot(ttr,range(psitrStoch0,2),'LineWidth',2)
legend('steady, \epsilon=0.01','stochastic, \kappa=10^{-6}','both', '\kappa=10^{-6}')
xlabel('time','fontsize',14)
set(gca,'fontsize',14)
ylabel('std(\Psi)','fontsize',14)
title('Standard deviation of Ensembles in \Psi, E=0.125, outer')

figure
subplot(3,1,1)
plot(ttr,psitrSteady1,'Color',[0.5 0.5 0.5])
hold on
plot(ttr,mean(psitrSteady1,2),'k','LineWidth',2)
plot(ttr,mean(psitrSteady1,2)+std(psitrSteady1,0,2),'k-.','LineWidth',2)
plot(ttr,mean(psitrSteady1,2)-std(psitrSteady1,0,2),'k-.','LineWidth',2)
ylim([0 1.8e-3])
set(gca,'fontsize',14)
title('E=0.125, Steady \epsilon=0.01','fontsize',14)
subplot(3,1,2)
plot(ttr,psitrBoth1,'Color',[0.5 0.5 0.5])
hold on
plot(ttr,mean(psitrBoth1,2),'k','LineWidth',2)
plot(ttr,mean(psitrBoth1,2)+std(psitrBoth1,0,2),'k-.','LineWidth',2)
plot(ttr,mean(psitrBoth1,2)-std(psitrBoth1,0,2),'k-.','LineWidth',2)
ylim([0 1.8e-3])
set(gca,'fontsize',14)
title('Steady and Stochastic','fontsize',14)
subplot(3,1,3)
plot(ttr,psitrStoch1,'Color',[0.5 0.5 0.5])
hold on
plot(ttr,mean(psitrStoch1,2),'k','LineWidth',2)
plot(ttr,mean(psitrStoch1,2)+std(psitrStoch1,0,2),'k-.','LineWidth',2)
plot(ttr,mean(psitrStoch1,2)-std(psitrStoch1,0,2),'k-.','LineWidth',2)
ylim([0 1.8e-3])
title('Stochastic, \kappa=10^{-6}','fontsize',14)
set(gca,'fontsize',14)
xlabel('time','fontsize',14)
%% E=0.02 outer, starting at x,y=0, z=0.1
load('E0.02sigma0.042eps0.01y0z0.1SphereV3.mat')
psitrBoth1=mypsi(R,E,xtr,ytr,ztr);
ztrBoth1=ztr;
rtrBoth1=r;
[phiBoth1,~]=cart2pol(ztrBoth1-0.5,rtrBoth1-0.5);
load('E0.02sigma0.042eps0y0z0.1SphereV3.mat')
psitrStoch1=mypsi(R,E,xtr,ytr,ztr);
ztrStoch1=ztr;
rtrStoch1=r;
[phiStoch1,~]=cart2pol(ztrStoch1-0.5,rtrStoch1-0.5);
load('E0.02sigma0eps0.01y0z0.1Sphere.mat')
psitrSteady1=mypsi(R,E,xtr,ytr,ztr);
ztrSteady1=ztr;
rtrSteady1=r;
[phiSteady1,~]=cart2pol(ztrSteady1-0.5,rtrSteady1-0.5);

% load('E0.02sigma0.042eps0.01y0.1z0.1Sphere.mat')
% psitrBoth1=mypsi(R,E,xtr,ytr,ztr);
% load('E0.02sigma0.042eps0y0.1z0.1Sphere.mat')
% psitrStoch1=mypsi(R,E,xtr,ytr,ztr);
% load('E0.02sigma0eps0.01y0.1z0.1Sphere.mat')
% psitrSteady1=mypsi(R,E,xtr,ytr,ztr);

%load('E0.125sigma0.013eps0r0.5Sphere.mat')
%psitrStochSmall=mypsi(R,E,xtr,ytr,ztr);

figure
plot(ttr,range(psitrSteady1.'),'LineWidth',2)
hold all
plot(ttr,range(psitrStoch1.'),'LineWidth',2)
plot(ttr,range(psitrBoth1.'),'LineWidth',2)
%plot(ttr,range(psitrStochSmall,2),'LineWidth',2)
legend('steady, \epsilon=0.01','stochastic, \kappa=10^{-6}','both', '\kappa=10^{-7}')
xlabel('time','fontsize',14)
set(gca,'fontsize',14)
ylabel('range in \Psi','fontsize',14)
title('Range of Ensembles in \Psi, E=0.02, centerline')


figure
subplot(3,1,1)
plot(ttr,psitrSteady1,'Color',[0.5 0.5 0.5])
hold on
plot(ttr,mean(psitrSteady1,2),'k','LineWidth',2)
plot(ttr,mean(psitrSteady1,2)+std(psitrSteady1,0,2),'k-.','LineWidth',2)
plot(ttr,mean(psitrSteady1,2)-std(psitrSteady1,0,2),'k-.','LineWidth',2)
ylim([0 1.8e-3])
set(gca,'fontsize',14)
title('E=0.02, Steady \epsilon=0.01','fontsize',14)
subplot(3,1,2)
plot(ttr,psitrBoth1,'Color',[0.5 0.5 0.5])
hold on
plot(ttr,mean(psitrBoth1,2),'k','LineWidth',2)
plot(ttr,mean(psitrBoth1,2)+std(psitrBoth1,0,2),'k-.','LineWidth',2)
plot(ttr,mean(psitrBoth1,2)-std(psitrBoth1,0,2),'k-.','LineWidth',2)
ylim([0 1.8e-3])
set(gca,'fontsize',14)
title('Steady and Stochastic','fontsize',14)
subplot(3,1,3)
plot(ttr,psitrStoch1,'Color',[0.5 0.5 0.5])
hold on
plot(ttr,mean(psitrStoch1,2),'k','LineWidth',2)
plot(ttr,mean(psitrStoch1,2)+std(psitrStoch1,0,2),'k-.','LineWidth',2)
plot(ttr,mean(psitrStoch1,2)-std(psitrStoch1,0,2),'k-.','LineWidth',2)
ylim([0 1.8e-3])
title('Stochastic, \kappa=10^{-6}','fontsize',14)
set(gca,'fontsize',14)
xlabel('time','fontsize',14)
%% wtf stochastic stuff
load('E0.125sigma0.042eps0r0.9Sphere.mat')
psitrStoch1=mypsi(R,E,xtr,ytr,ztr);
xtr1=xtr; ytr1=ytr; ztr1=ztr;
load('E0.125sigma0.042eps0y0.05z0.1SphereV2.mat')
psitrStoch2=mypsi(R,E,xtr,ytr,ztr);
xtr2=xtr; ytr2=ytr; ztr2=ztr;

figure; plot(mean(xtr1,2)); hold all; plot(mean(xtr2,2)); title('mean x')
figure; plot(range(xtr1,2)); hold all; plot(range(xtr2,2)); title('range x')
figure; plot(std(xtr1,0,2)); hold all; plot(std(xtr2,0,2)); title('std x')

figure; plot(mean(psitrStoch1,2)-mean(psitrStoch1(1,:))); hold all; plot(mean(psitrStoch2,2)-mean(psitrStoch2(1,:))); title('change in mean \psi from t=0')
figure; plot(range(psitrStoch1,2)-range(psitrStoch1(1,:))); hold all; plot(range(psitrStoch2,2)-range(psitrStoch2(1,:))); title('change in range \psi from t=0')
figure; plot(std(psitrStoch1,0,2)-std(psitrStoch1(1,:),0,2)); hold all; plot(std(psitrStoch2,0,2)-std(psitrStoch2(1,:),0,2)); title('change in std \psi from t=0')

%% full set load, figs
pert=4; x0=-0.9;%-0.5;
Eks=[0.25 0.125 0.02 0.0005];

load 'initialSpheres.mat'

subi=0;
figure;
for Eki=2:3%:2
    E=Eks(Eki)
    subi=subi+1
    for rInit=[0.1 0.4]%0:0.1:0.5%[0.1 0.3 0.7 0.9]%[0.5]% 0.9]
        if rInit==0.4
            subi=subi+1
        end
        subplot(2,2,subi); hold on;
        ni=1;
        for pertType=1:6%[1 3]
            switch pertType
                case 5
                    sigma=0.042;
                    eps=0.0;
                    kappa=1e-6;
                    ftitle1s=strcat('E',num2str(E),'sigma',num2str(sigma),'eps',num2str(eps),'r',num2str(rInit),'S.mat');
                    c1=[0 0 1];
                case 1
                    sigma=0.0;
                    eps=0.01;
                    kappa=0;
                    ftitle1s=strcat('E',num2str(E),'sigma',num2str(sigma),'eps',num2str(eps),'r',num2str(rInit),'S.mat');%Sb
                    c1=[1 0 0];
                case 3
                    sigma=0.013;%42;
                    eps=0.01;
                    kappa=1e-7;%6;
                    ftitle1s=strcat('E',num2str(E),'sigma',num2str(sigma),'eps',num2str(eps),'r',num2str(rInit),'S.mat');%Sb
                    c1=[0 0.5 0];
                case 4
                    sigma=0.013;
                    eps=0.0;
                    kappa=1e-7;
                    ftitle1s=strcat('E',num2str(E),'sigma',num2str(sigma),'eps',num2str(eps),'r',num2str(rInit),'S.mat');
                    c1=[0 1 1];
                case 6
                    sigma=0.13;
                    eps=0.0;
                    kappa=1e-5;
                    ftitle1s=strcat('E',num2str(E),'sigma',num2str(sigma),'eps',num2str(eps),'r',num2str(rInit),'S.mat');
                    c1=[0.1 0.2 0.5];
                case 2
                    sigma=0;
                    eps=0.08;
                    kappa=0;
                    ftitle1s=strcat('E',num2str(E),'sigma',num2str(sigma),'eps',num2str(eps),'r',num2str(rInit),'S.mat');%Sb
                    c1=[1 0 1];
            end
    
            %ftitle1s=strcat('E',num2str(E),'sigma',num2str(sigma),'eps',num2str(eps),'r',num2str(rInit),'Sb.mat');
            load(ftitle1s)
            %trajectorySlice
            psitr=mypsi(R,E,xtr,ytr,ztr);
           plot(ttr,range(psitr.'),'linewidth',2,'Color',c1); hold all
           %plot(ttr,std(psitr.')); hold all
           legend1(ni)=cellstr(strcat('\epsilon=',num2str(eps),', \kappa=',num2str(kappa)));
            ni=ni+1;
            
        end %for pert
        %title(strcat('Range in \psi, E=',num2str(E),', r=',num2str(rInit)),'fontsize',18)
        if subi==1
        xlabel('(a)','fontsize',16)
        set(gca,'XTickLabel',{})
        ylabel('range in \psi','fontsize',16)
        elseif subi==2
            xlabel('(b)','fontsize',16)
            set(gca,'XTickLabel',{})
            set(gca,'YTickLabel',{})
        elseif subi==3
            xlabel({'time';'(c)'},'fontsize',16)
            ylabel('range in \psi','fontsize',16)
        elseif subi==4
            xlabel({'time';'(d)'},'fontsize',16)
            set(gca,'YTickLabel',{})
            legend(legend1)
        end
        
        set(gca,'fontsize',16)
        axis([0 10000 0 2.4e-3])
    end %for rInit

end %for E