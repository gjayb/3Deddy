%was figures for OSM poster (2014); now using for some paper figures (2018)
%%%%%
%% rc paper figure 1
%parameters and grid
E=0.125;
R=1;
x=-1:0.005:1;
z=0:0.001:1;
[xx,zz]=meshgrid(x,z);
%constants
Eratio=sqrt(E);%ratio of Ekman layer thickness to dimensional tank depth
s=sin(0.5/Eratio);
c=cos(0.5/Eratio);
s2=sinh(0.5/Eratio);
c2=cosh(0.5/Eratio);
A=-0.5*c*s2/((s^2)*(c2^2)+(c^2)*(s2^2));
B=0.5*c2*s/((s^2)*(c2^2)+(c^2)*(s2^2));
D=A*(s*s2-c*c2)+B*(c*c2+s*s2);
%coordinates
r=abs(xx);
zeta=(zz-0.5)/Eratio;
%streamfunction; 
Rfun=0.5.*(r.^2).*((R-r).^2);
Ffun=(A.*(sin(zeta).*sinh(zeta)-cos(zeta).*cosh(zeta))...
    +B.*(cos(zeta).*cosh(zeta)+sin(zeta).*sinh(zeta))-D);
psi=-Eratio.*Rfun.*Ffun;

 figure
 subplot(2,3,1)
 contourf(xx,zz,psi,0:max(psi(:))/10:max(psi(:)));%E=0.25, use 0:0.00008:max(max(psi))
 colormap(gray)
 caxis([-0.0003 max(max(psi))])
 %axis equal
 axis tight
 ylabel('Z','fontsize',14)
  set(gca,'ytick',[0 0.5 1])
 title(strcat('E=',num2str(E)),'fontsize',16)
 xlabel({'Y','(a)'},'fontsize',14)
 set(gca,'xtick',[-1 0 1])
 set(gca,'fontsize',12)
freezeColors
hold on; plot([-0.5 0.5],[0.5 0.5],'b.')

E=0.02;
Eratio=sqrt(E);%ratio of Ekman layer thickness to dimensional tank depth
s=sin(0.5/Eratio);
c=cos(0.5/Eratio);
s2=sinh(0.5/Eratio);
c2=cosh(0.5/Eratio);
A=-0.5*c*s2/((s^2)*(c2^2)+(c^2)*(s2^2));
B=0.5*c2*s/((s^2)*(c2^2)+(c^2)*(s2^2));
D=A*(s*s2-c*c2)+B*(c*c2+s*s2);
%coordinates
r=abs(xx);
zeta=(zz-0.5)/Eratio;
%streamfunction; 
Rfun=0.5.*(r.^2).*((R-r).^2);
Ffun=(A.*(sin(zeta).*sinh(zeta)-cos(zeta).*cosh(zeta))...
    +B.*(cos(zeta).*cosh(zeta)+sin(zeta).*sinh(zeta))-D);
psi=-Eratio.*Rfun.*Ffun;

subplot(2,3,2)
contourf(xx,zz,psi,0:max(psi(:))/10:max(psi(:)));%E=0.01, use [16.1 psi(501,301)*1e4 15.2  14.2 13 11.5 10 8 6 4 1 0].*1e-4
colormap(gray)
caxis([-0.0004 max(max(psi))])
 %axis equal
 axis tight
title(strcat('E=',num2str(E)),'fontsize',16)
 set(gca,'ytick',[])
  xlabel({'Y','(b)'},'fontsize',14)
 set(gca,'xtick',[-1 0 1])
 set(gca,'fontsize',12)
freezeColors
hold on; plot([-0.5 0.5],[0.5 0.5],'b.')
 
 E=0.0005;
Eratio=sqrt(E);%ratio of Ekman layer thickness to dimensional tank depth
s=sin(0.5/Eratio);
c=cos(0.5/Eratio);
s2=sinh(0.5/Eratio);
c2=cosh(0.5/Eratio);
A=-0.5*c*s2/((s^2)*(c2^2)+(c^2)*(s2^2));
B=0.5*c2*s/((s^2)*(c2^2)+(c^2)*(s2^2));
D=A*(s*s2-c*c2)+B*(c*c2+s*s2);
%coordinates
r=abs(xx);
zeta=(zz-0.5)/Eratio;
%streamfunction; 
Rfun=0.5.*(r.^2).*((R-r).^2);
Ffun=(A.*(sin(zeta).*sinh(zeta)-cos(zeta).*cosh(zeta))...
    +B.*(cos(zeta).*cosh(zeta)+sin(zeta).*sinh(zeta))-D);
psi=-Eratio.*Rfun.*Ffun;

subplot(2,3,3)
contourf(xx,zz,psi,0:max(psi(:))/10:max(psi(:)));%E=0.001, use [5.05 4.98 psi(501,301)*1e4 4.91 4.6 4.1 3.3 2.5 1.7 0.9 0.1 0].*1e-4
colormap(gray)
caxis([-0.00014 max(max(psi))])
 %axis equal
 axis tight
title(strcat('E=',num2str(E)),'fontsize',16)
 set(gca,'ytick',[])
  xlabel({'Y';'(c)'},'fontsize',14)
 set(gca,'xtick',[-1 0 1])
 set(gca,'fontsize',12)
freezeColors
hold on; plot(-[0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5],[0.93 0.86 0.789 0.719 0.649 0.579 0.5 0.421 0.351 0.281 0.211 0.14 0.07],'b.')
plot([0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5],[0.93 0.86 0.789 0.719 0.649 0.579 0.5 0.421 0.351 0.281 0.211 0.14 0.07],'b.')
 
 %azimuthal perturbation
yin=-1:0.01:1;
xin=-1:0.01:1;
[xx,yy,zz]=meshgrid(xin,yin,1);
R=1;
[theta,r]=cart2pol(xx,yy);
E=0.25;
x0=-0.5;
Eratio=sqrt(E);
r2=sqrt((xx-x0).^2 +yy.^2);
psi1=-sinh(zz/Eratio).*(R.^2 -r.^2).*(4*R.^2 -r2.^2)./sinh(1/Eratio);
psi1(r>R)=NaN;

subplot(2,3,5)
contourf(squeeze(xx),squeeze(yy),squeeze(psi1),20)
colormap(darkb2r(-8,1))
axis equal
axis tight
set(gca,'color',[0.5 0.5 0.5])
 xlabel({'X';'(d)'},'fontsize',14)
 set(gca,'xtick',[-1 0 1])
 ylabel('Y','fontsize',14)
 set(gca,'ytick',[-1 0 1])
 set(gca,'fontsize',12)

 
%%
%%1: background overturning streamfunction; E=1, E=0.01, E=0.0001

%parameters and grid
E=0.125;
R=1;
x=-1:0.005:1;
z=0:0.001:1;
[xx,zz]=meshgrid(x,z);
%constants
Eratio=sqrt(E);%ratio of Ekman layer thickness to dimensional tank depth
s=sin(0.5/Eratio);
c=cos(0.5/Eratio);
s2=sinh(0.5/Eratio);
c2=cosh(0.5/Eratio);
A=-0.5*c*s2/((s^2)*(c2^2)+(c^2)*(s2^2));
B=0.5*c2*s/((s^2)*(c2^2)+(c^2)*(s2^2));
D=A*(s*s2-c*c2)+B*(c*c2+s*s2);
%coordinates
r=abs(xx);
zeta=(zz-0.5)/Eratio;
%streamfunction; 
Rfun=0.5.*(r.^2).*((R-r).^2);
Ffun=(A.*(sin(zeta).*sinh(zeta)-cos(zeta).*cosh(zeta))...
    +B.*(cos(zeta).*cosh(zeta)+sin(zeta).*sinh(zeta))-D);
psi=-Eratio.*Rfun.*Ffun;

 figure
 subplot('position',[0.06 0.1 0.29 0.8])
 contourf(xx,zz,psi,0:max(psi(:))/10:max(psi(:)));%E=0.25, use 0:0.00008:max(max(psi))
 colormap(gray)
 caxis([-0.0003 max(max(psi))])
 %axis equal
 axis tight
 xlabel('Y','fontsize',14)
 ylabel('Z','fontsize',14)
 title(strcat('E=',num2str(E)),'fontsize',16)
 set(gca,'xtick',[-1 0 1])
 set(gca,'ytick',[0 0.5 1])

E=0.02;
Eratio=sqrt(E);%ratio of Ekman layer thickness to dimensional tank depth
s=sin(0.5/Eratio);
c=cos(0.5/Eratio);
s2=sinh(0.5/Eratio);
c2=cosh(0.5/Eratio);
A=-0.5*c*s2/((s^2)*(c2^2)+(c^2)*(s2^2));
B=0.5*c2*s/((s^2)*(c2^2)+(c^2)*(s2^2));
D=A*(s*s2-c*c2)+B*(c*c2+s*s2);
%coordinates
r=abs(xx);
zeta=(zz-0.5)/Eratio;
%streamfunction; 
Rfun=0.5.*(r.^2).*((R-r).^2);
Ffun=(A.*(sin(zeta).*sinh(zeta)-cos(zeta).*cosh(zeta))...
    +B.*(cos(zeta).*cosh(zeta)+sin(zeta).*sinh(zeta))-D);
psi=-Eratio.*Rfun.*Ffun;

subplot('position',[0.38 0.1 0.29 0.8])
contourf(xx,zz,psi,0:max(psi(:))/10:max(psi(:)));%E=0.01, use [16.1 psi(501,301)*1e4 15.2  14.2 13 11.5 10 8 6 4 1 0].*1e-4
colormap(gray)
caxis([-0.0004 max(max(psi))])
 %axis equal
 axis tight
title(strcat('E=',num2str(E)),'fontsize',16)
 set(gca,'xtick',[])
 set(gca,'ytick',[])
 
 E=0.0005;
Eratio=sqrt(E);%ratio of Ekman layer thickness to dimensional tank depth
s=sin(0.5/Eratio);
c=cos(0.5/Eratio);
s2=sinh(0.5/Eratio);
c2=cosh(0.5/Eratio);
A=-0.5*c*s2/((s^2)*(c2^2)+(c^2)*(s2^2));
B=0.5*c2*s/((s^2)*(c2^2)+(c^2)*(s2^2));
D=A*(s*s2-c*c2)+B*(c*c2+s*s2);
%coordinates
r=abs(xx);
zeta=(zz-0.5)/Eratio;
%streamfunction; 
Rfun=0.5.*(r.^2).*((R-r).^2);
Ffun=(A.*(sin(zeta).*sinh(zeta)-cos(zeta).*cosh(zeta))...
    +B.*(cos(zeta).*cosh(zeta)+sin(zeta).*sinh(zeta))-D);
psi=-Eratio.*Rfun.*Ffun;

subplot('position',[0.7 0.1 0.29 0.8])
contourf(xx,zz,psi,0:max(psi(:))/10:max(psi(:)));%E=0.001, use [5.05 4.98 psi(501,301)*1e4 4.91 4.6 4.1 3.3 2.5 1.7 0.9 0.1 0].*1e-4
colormap(gray)
caxis([-0.00014 max(max(psi))])
 %axis equal
 axis tight
title(strcat('E=',num2str(E)),'fontsize',16)
 set(gca,'xtick',[])
 set(gca,'ytick',[])
 
 clear
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%
 
%2: Poincare slices to compare to Navier-Stokes sim ones
figure
%subplot('position',[0.06 0.1 0.2 0.8])
load('E0.25P4e0.01x0-0.5.mat')
hold all;
for itr=1:size(ztr,2);%[15 59 70 85];%:size(xtr,2)
             jj=find(xtr(1:end-1,itr).*xtr(2:end,itr)<0); % look for 
%points along the trajectory where x changes sign
             plot(ytr(jj,itr),ztr(jj,itr),'k.','MarkerSize',1);
end
xlabel('Y','fontsize',14)
 ylabel('Z','fontsize',14)
 title(strcat('E=',num2str(E)),'fontsize',16)
 set(gca,'xtick',[-1 0 1])
 set(gca,'ytick',[0 0.5 1])
 axis equal
axis tight

figure
hold all
for itr=1:size(ztr,2);%[15 59 70 85];%:size(xtr,2)
             jj=find(theta(1:end-1,itr)<(pi/2)+0.0075 & theta(1:end-1,itr)>(pi/2)-0.0075); % look for 
%points where we are near the y-axis
            jj2=find(theta(1:end-1,itr)<(-pi/2)+0.0075 & theta(1:end-1,itr)>(-pi/2)-0.0075);
             plot(ytr(jj,itr),ztr(jj,itr),'k.','MarkerSize',1);
             plot(ytr(jj2,itr),ztr(jj2,itr),'k.','MarkerSize',1);
end
xlabel('Y','fontsize',14)
 ylabel('Z','fontsize',14)
 title(strcat('E=',num2str(E)),'fontsize',16)
 set(gca,'xtick',[-1 0 1])
 set(gca,'ytick',[0 0.5 1])
 axis equal
axis tight
 
figure
hold all
for itr=1:size(ztr,2);%[15 59 70 85];%:size(xtr,2)
             jj=find(xtr(1:end-1,itr)<0.0075 & xtr(1:end-1,itr)>-0.0075); % look for 
%points where we are near the y-axis
            jj2=find(xtr(1:end-1,itr)<+0.0075 & xtr(1:end-1,itr)>-0.0075);
             plot(ytr(jj,itr),ztr(jj,itr),'k.','MarkerSize',1);
             plot(ytr(jj2,itr),ztr(jj2,itr),'k.','MarkerSize',1);
end
xlabel('Y','fontsize',14)
 ylabel('Z','fontsize',14)
 title(strcat('E=',num2str(E)),'fontsize',16)
 set(gca,'xtick',[-1 0 1])
 set(gca,'ytick',[0 0.5 1])
 axis equal
axis tight

% subplot('position',[0.3 0.1 0.2 0.8])
figure
load('E0.125P4e0.05x0-0.9.mat')
hold all;
for itr=1:size(ztr,2);%[15 59 70 85];%:size(xtr,2)
             jj=find(xtr(1:end-1,itr).*xtr(2:end,itr)<0); % look for 
%points along the trajectory where x changes sign
             plot(ytr(jj,itr),ztr(jj,itr),'k.','MarkerSize',1);

end
 title(strcat('E=',num2str(E)),'fontsize',16)
 set(gca,'xtick',[])
 set(gca,'ytick',[])
 axis equal
axis tight
 
figure
hold all
for itr=1:size(ztr,2);%[15 59 70 85];%:size(xtr,2)
             jj=find(theta(1:end-1,itr)<(pi/2)+0.0075 & theta(1:end-1,itr)>(pi/2)-0.0075); % look for 
%points where we are near the y-axis
            jj2=find(theta(1:end-1,itr)<(-pi/2)+0.0075 & theta(1:end-1,itr)>(-pi/2)-0.0075);
             plot(ytr(jj,itr),ztr(jj,itr),'k.','MarkerSize',1);
             plot(ytr(jj2,itr),ztr(jj2,itr),'k.','MarkerSize',1);
end
xlabel('Y','fontsize',14)
 ylabel('Z','fontsize',14)
 title(strcat('E=',num2str(E)),'fontsize',16)
 set(gca,'xtick',[-1 0 1])
 set(gca,'ytick',[0 0.5 1])
 axis equal
axis tight
 
figure
hold all
for itr=1:size(ztr,2);%[15 59 70 85];%:size(xtr,2)
             jj=find(xtr(1:end-1,itr)<0.0075 & xtr(1:end-1,itr)>-0.0075); % look for 
%points where we are near the y-axis
            jj2=find(xtr(1:end-1,itr)<+0.0075 & xtr(1:end-1,itr)>-0.0075);
             plot(ytr(jj,itr),ztr(jj,itr),'k.','MarkerSize',1);
             plot(ytr(jj2,itr),ztr(jj2,itr),'k.','MarkerSize',1);
end
xlabel('Y','fontsize',14)
 ylabel('Z','fontsize',14)
 title(strcat('E=',num2str(E)),'fontsize',16)
 set(gca,'xtick',[-1 0 1])
 set(gca,'ytick',[0 0.5 1])
 axis equal
axis tight

% subplot('position',[0.54 0.1 0.2 0.8])
figure
load('E0.02P4e0.1x0-0.5.mat')
hold all;
for itr=1:size(ztr,2);%[15 59 70 85];%:size(xtr,2)
             jj=find(xtr(1:end-1,itr).*xtr(2:end,itr)<0); % look for 
%points along the trajectory where x changes sign
             plot(ytr(jj,itr),ztr(jj,itr),'k.','MarkerSize',1);
end
 title(strcat('E=',num2str(E)),'fontsize',16)
 set(gca,'xtick',[])
 set(gca,'ytick',[])
 axis equal
axis tight
 
figure
hold all
for itr=1:size(ztr,2);%[15 59 70 85];%:size(xtr,2)
             jj=find(theta(1:end-1,itr)<(pi/2)+0.0075 & theta(1:end-1,itr)>(pi/2)-0.0075); % look for 
%points where we are near the y-axis
            jj2=find(theta(1:end-1,itr)<(-pi/2)+0.0075 & theta(1:end-1,itr)>(-pi/2)-0.0075);
             plot(ytr(jj,itr),ztr(jj,itr),'k.','MarkerSize',1);
             plot(ytr(jj2,itr),ztr(jj2,itr),'k.','MarkerSize',1);
end
xlabel('Y','fontsize',14)
 ylabel('Z','fontsize',14)
 title(strcat('E=',num2str(E)),'fontsize',16)
 set(gca,'xtick',[-1 0 1])
 set(gca,'ytick',[0 0.5 1])
 axis equal
axis tight
 
figure
hold all
for itr=1:size(ztr,2);%[15 59 70 85];%:size(xtr,2)
             jj=find(xtr(1:end-1,itr)<0.0075 & xtr(1:end-1,itr)>-0.0075); % look for 
%points where we are near the y-axis
            jj2=find(xtr(1:end-1,itr)<+0.0075 & xtr(1:end-1,itr)>-0.0075);
             plot(ytr(jj,itr),ztr(jj,itr),'k.','MarkerSize',1);
             plot(ytr(jj2,itr),ztr(jj2,itr),'k.','MarkerSize',1);
end
xlabel('Y','fontsize',14)
 ylabel('Z','fontsize',14)
 title(strcat('E=',num2str(E)),'fontsize',16)
 set(gca,'xtick',[-1 0 1])
 set(gca,'ytick',[0 0.5 1])
 axis equal
axis tight

% subplot('position',[0.78 0.1 0.2 0.8])
figure
load('E0.0005P4e0.01x0-0.9.mat')
hold all;
for itr=1:size(ztr,2);%[15 59 70 85];%:size(xtr,2)
             jj=find(xtr(1:end-1,itr).*xtr(2:end,itr)<0); % look for 
%points along the trajectory where x changes sign
             plot(ytr(jj,itr),ztr(jj,itr),'k.','MarkerSize',1);

end
 title(strcat('E=',num2str(E)),'fontsize',16)
 set(gca,'xtick',[])
 set(gca,'ytick',[])
 axis equal
axis tight

figure
hold all
for itr=1:size(ztr,2);%[15 59 70 85];%:size(xtr,2)
             jj=find(theta(1:end-1,itr)<(pi/2)+0.0075 & theta(1:end-1,itr)>(pi/2)-0.0075); % look for 
%points where we are near the y-axis
            jj2=find(theta(1:end-1,itr)<(-pi/2)+0.0075 & theta(1:end-1,itr)>(-pi/2)-0.0075);
             plot(ytr(jj,itr),ztr(jj,itr),'k.','MarkerSize',1);
             plot(ytr(jj2,itr),ztr(jj2,itr),'k.','MarkerSize',1);
end
xlabel('Y','fontsize',14)
 ylabel('Z','fontsize',14)
 title(strcat('E=',num2str(E)),'fontsize',16)
 set(gca,'xtick',[-1 0 1])
 set(gca,'ytick',[0 0.5 1])
 axis equal
axis tight
 
figure
hold all
for itr=1:size(ztr,2);%[15 59 70 85];%:size(xtr,2)
             jj=find(xtr(1:end-1,itr)<0.0075 & xtr(1:end-1,itr)>-0.0075); % look for 
%points where we are near the y-axis
            jj2=find(xtr(1:end-1,itr)<+0.0075 & xtr(1:end-1,itr)>-0.0075);
             plot(ytr(jj,itr),ztr(jj,itr),'k.','MarkerSize',1);
             plot(ytr(jj2,itr),ztr(jj2,itr),'k.','MarkerSize',1);
end
xlabel('Y','fontsize',14)
 ylabel('Z','fontsize',14)
 title(strcat('E=',num2str(E)),'fontsize',16)
 set(gca,'xtick',[-1 0 1])
 set(gca,'ytick',[0 0.5 1])
 axis equal
axis tight

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

%3: Trajectories to demonstrate effects of noise and steady perturbation
figure
load('E0125eps0sigma0fig.mat')
subplot('position',[0.06 0.1 0.2 0.8])
plot3(xtr(1:40000,1:4),ytr(1:40000,1:4),ztr(1:40000,1:4))
axis([-1 1 -1 1 0 1])
set(gca,'ztick',[0 0.5 1])
set(gca,'xtick',[-1 0 1])
set(gca,'ytick',[-1 0 1])
zlabel('z','fontsize',14)
grid on
title({'Background Case:','symmetric, deterministic'},'fontsize',14)

load('E0125eps005sigma0fig.mat')
subplot('position',[0.3 0.1 0.2 0.8])
plot3(xtr(1:50000,1:4),ytr(1:50000,1:4),ztr(1:50000,1:4))
axis([-1 1 -1 1 0 1])
set(gca,'ztick',[0 0.5 1])
set(gca,'xtick',[-1 0 1])
set(gca,'ytick',[-1 0 1])
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
set(gca,'zticklabel',[])
grid on
title({'Background plus steady','symmetry-breaking perturbation'},'fontsize',16)

load('E0125eps0sigma005fig.mat')
subplot('position',[0.54 0.1 0.2 0.8])
plot3(xtr(1:40000,1:4),ytr(1:40000,1:4),ztr(1:40000,1:4))
axis([-1 1 -1 1 0 1])
set(gca,'ztick',[0 0.5 1])
set(gca,'xtick',[-1 0 1])
set(gca,'ytick',[-1 0 1])
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
set(gca,'zticklabel',[])
grid on
title({'Background plus','random-walk noise'},'fontsize',16)


 subplot('position',[0.78 0.1 0.2 0.8])
 load('E0125eps005sigma005fig.mat')
plot3(xtr(1:30000,1:4),ytr(1:30000,1:4),ztr(1:30000,1:4))
axis([-1 1 -1 1 0 1])
set(gca,'ztick',[0 0.5 1])
set(gca,'xtick',[-1 0 1])
set(gca,'ytick',[-1 0 1])
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
set(gca,'zticklabel',[])
grid on
title({'Background plus noise','and Steady Perturbation'},'fontsize',16)
%%%%%%%%%%%%%%%%%%%%%%%
%%

%4: Lyapunov exponents to show regular vs chaotic regions.
figure
load('ftleE025.mat')
i=4;
[~,Fp,Fn,Fm,~]=ftleCalc(xin3b,yin3b,zin3b,xtr(i+1,:),ytr(i+1,:),ztr(i+1,:),i);
subplot('position',[0.16 0.1 0.1 0.8])
pcolor(yin3b(2:end-1,1,2:end-1),zin3b(2:end-1,1,2:end-1),squeeze(Fp))
shading flat
axis equal
axis tight
colormap(spring)
colorbar('SouthOutside')
title('E=0.25','fontsize',16)
xlabel('y','fontsize',14)
ylabel('z','fontsize',14)
set(gca,'xtick',[0.25 0.75])
set(gca,'ytick',[0.25 0.75])

subplot('position',[0.40 0.1 0.1 0.8])
load('ftleP4E0125.mat')
i=4;
[~,Fp,Fn,Fm,~]=ftleCalc(xin3b,yin3b,zin3b,xtr(i+1,:),ytr(i+1,:),ztr(i+1,:),i);
pcolor(yin3b(2:end-1,1,2:end-1),zin3b(2:end-1,1,2:end-1),squeeze(Fp))
shading flat
axis equal
axis tight
colormap(spring)
colorbar('SouthOutside')
title('E=0.125','fontsize',16)
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])

subplot('position',[0.64 0.1 0.1 0.8])
load('ftleP4E002.mat')
i=5;
[~,Fp,Fn,Fm,~]=ftleCalc(xin3b,yin3b,zin3b,xtr(i+1,:),ytr(i+1,:),ztr(i+1,:),i);
pcolor(yin3b(2:end-1,1,2:end-1),zin3b(2:end-1,1,2:end-1),squeeze(Fp))
shading flat
axis equal
axis tight
colormap(spring)
colorbar('SouthOutside')
title('E=0.02','fontsize',16)
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])

subplot('position',[0.88 0.1 0.1 0.8])
load('ftleP4E00005.mat')
i=5;
[~,Fp,Fn,Fm,~]=ftleCalc(xin3b,yin3b,zin3b,xtr(i+1,:),ytr(i+1,:),ztr(i+1,:),i);
pcolor(yin3b(2:end-1,1,2:end-1),zin3b(2:end-1,1,2:end-1),squeeze(Fp))
shading flat
axis equal
axis tight
colormap(spring)
colorbar('SouthOutside')
title('E=0.0005','fontsize',16)
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])

%%%%%%%%%%%%%%%%%%%%%%
%%

%5: Lyapunov exponents 2: show more accurate (shorter time) values

figure
load('ftleShortP4E025.mat')
i=3;
[~,Fp,Fn,Fm,~]=ftleCalc(xin3b,yin3b,zin3b,xtr(i+1,:),ytr(i+1,:),ztr(i+1,:),i);
subplot('position',[0.06 0.1 0.2 0.8])
pcolor(yin3b(2:end-1,1,2:end-1),zin3b(2:end-1,1,2:end-1),squeeze(Fn))
shading flat
axis equal
axis tight
colormap(cool)
colorbar('SouthOutside')
caxis([-3 0])
title('E=0.25','fontsize',16)
xlabel('y','fontsize',14)
ylabel('z','fontsize',14)
set(gca,'xtick',[0.25 0.75])
set(gca,'ytick',[0.25 0.75])

subplot('position',[0.3 0.1 0.2 0.8])
load('ftleShortP4E0125.mat')
i=3;
[~,Fp,Fn,Fm,~]=ftleCalc(xin3b,yin3b,zin3b,xtr(i+1,:),ytr(i+1,:),ztr(i+1,:),i);
pcolor(yin3b(2:end-1,1,2:end-1),zin3b(2:end-1,1,2:end-1),squeeze(Fn))
shading flat
axis equal
axis tight
colormap(cool)
colorbar('SouthOutside')
caxis([-3 0])
title('E=0.125','fontsize',16)
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])

subplot('position',[0.54 0.1 0.2 0.8])
load('ftleShortP4E002.mat')
i=3;
[~,Fp,Fn,Fm,~]=ftleCalc(xin3b,yin3b,zin3b,xtr(i+1,:),ytr(i+1,:),ztr(i+1,:),i);
pcolor(yin3b(2:end-1,1,2:end-1),zin3b(2:end-1,1,2:end-1),squeeze(Fn))
shading flat
axis equal
axis tight
colormap(cool)
colorbar('SouthOutside')
caxis([-3 0])
title('E=0.02','fontsize',16)
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])

subplot('position',[0.78 0.1 0.2 0.8])
load('ftleShortP4E00005.mat')
i=3;
[~,Fp,Fn,Fm,~]=ftleCalc(xin3b,yin3b,zin3b,xtr(i+1,:),ytr(i+1,:),ztr(i+1,:),i);
pcolor(yin3b(2:end-1,1,2:end-1),zin3b(2:end-1,1,2:end-1),squeeze(Fn))
shading flat
axis equal
axis tight
colormap(cool)
colorbar('SouthOutside')
caxis([-3 0])
title('E=0.0005','fontsize',16)
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])

%%%%%%%%%%%%%%%%%%%%%
%%

%6: 2 types of filamentation geometry
load('sphereE0125P4e005Xi09Xc0Yc0.2Zc0.2.mat')
figure
subplot('position',[0.06 0.1 0.2 0.8])
i=20;
surf(reshape(xtr(i,:),size(xi)),reshape(ytr(i,:),size(xi)),reshape(ztr(i,:),size(xi)),C)
shading flat
axis equal
view(-81,-2)
title({'Cigar shape, t=100','view 1'},'fontsize',16)
ylabel('y','fontsize',14)
zlabel('z','fontsize',14)

subplot('position',[0.3 0.1 0.2 0.8])
i=20;
surf(reshape(xtr(i,:),size(xi)),reshape(ytr(i,:),size(xi)),reshape(ztr(i,:),size(xi)),C)
shading flat
axis equal
title('Cigar shape, view 2','fontsize',16)
ylabel('y','fontsize',14)
view(80,52)

load('sphereE0125P4e005Xi09Xc0Yc0.4Zc0.2.mat')
subplot('position',[0.54 0.1 0.2 0.8])
i=20;
surf(reshape(xtr(i,:),size(xi)),reshape(ytr(i,:),size(xi)),reshape(ztr(i,:),size(xi)),C)
shading flat
axis equal
title({'Pancake shape, t=100','view 1'},'fontsize',16)
ylabel('y','fontsize',14)
view(-71,16)

subplot('position',[0.78 0.1 0.2 0.8])
i=20;
surf(reshape(xtr(i,:),size(xi)),reshape(ytr(i,:),size(xi)),reshape(ztr(i,:),size(xi)),C)
shading flat
axis equal
title('Pancake shape, view 2','fontsize',16)
ylabel('y','fontsize',14)
view(34,20)
%%
%examining the azimuthal velocity including perturbation

yin=-1:0.01:1;
zin=0:0.01:1;
[xx,yy,zz]=meshgrid(0,yin,zin);
R=1;

E=0.25;
x0=-0.5;
eps=0.1;
Eratio=sqrt(E);%ratio of Ekman layer thickness to dimensional tank depth
s=sin(0.5/Eratio);
c=cos(0.5/Eratio);
s2=sinh(0.5/Eratio);
c2=cosh(0.5/Eratio);
A=-0.5*c*s2/((s^2)*(c2^2)+(c^2)*(s2^2));
B=0.5*c2*s/((s^2)*(c2^2)+(c^2)*(s2^2));
D=A*(s*s2-c*c2)+B*(c*c2+s*s2);

[theta, r]=cart2pol(xx,yy);
r2=sqrt((yy.^2)+(xx-x0).^2);
zeta=(zz-0.5)/Eratio;

ur=r.*((R-r).^2).*(A.*sin(zeta).*cosh(zeta)+B.*cos(zeta).*sinh(zeta));
vtheta=r.*((R-r).^2).*(0.5+B.*sin(zeta).*cosh(zeta)-A.*cos(zeta).*sinh(zeta));
u=cos(theta).*ur-r.*sin(theta).*vtheta;
v=sin(theta).*ur+r.*cos(theta).*vtheta;
vthetaB=(v.*cos(theta)-u.*sin(theta))./r;
% figure
% pcolor(yy,zz,squeeze(vtheta))
% colormap(darkb2r(-1,1))
% shading 'flat'
% colorbar
% axis equal
% axis tight
% figure
% pcolor(yy,zz,squeeze(vthetaB))
% colormap(darkb2r(-1,1))
% shading 'flat'
% colorbar
% axis equal
% axis tight
 uT=5.*u+eps.*(2.*yy.*(R - r.^2) + 2.*yy.*(4*R^2 - r2.^2)).*(sinh(zz./Eratio)./sinh(1./Eratio));
 vT=5.*v+eps.*(-2.*(-x0 + xx).*(R - r.^2) - 2.*xx.*(4*R^2 - r2.^2)).*(sinh(zz./Eratio)./sinh(1./Eratio));
vtheta2=(vT.*cos(theta)-uT.*sin(theta))./r;
vtheta3=5.*vtheta-2*eps.*sinh(zz./Eratio).*((R-r.^2)+(4*R-r2.^2)-x0*cos(theta).*(R-r.^2)./r)./sinh(1/Eratio);

figure
pcolor(yy,zz,squeeze(vtheta2))
colormap(darkb2r(-1,1))
shading 'flat'
colorbar
axis equal
axis tight
title('E 0.25, eps 0.1, numerical')
figure
pcolor(yy,zz,squeeze(vtheta3))
colormap(darkb2r(-1,1))
shading 'flat'
colorbar
axis equal
axis tight
title('E 0.25, eps 0.1, analytic')

E=0.25;
x0=-0.5;
eps=0.05;

Eratio=sqrt(E);
s=sin(0.5/Eratio);
c=cos(0.5/Eratio);
s2=sinh(0.5/Eratio);
c2=cosh(0.5/Eratio);
A=-0.5*c*s2/((s^2)*(c2^2)+(c^2)*(s2^2));
B=0.5*c2*s/((s^2)*(c2^2)+(c^2)*(s2^2));
D=A*(s*s2-c*c2)+B*(c*c2+s*s2);
r2=sqrt((yy.^2)+(xx-x0).^2);
zeta=(zz-0.5)/Eratio;
vtheta=r.*((R-r).^2).*(0.5+B.*sin(zeta).*cosh(zeta)-A.*cos(zeta).*sinh(zeta));
vtheta4=5.*vtheta-2*eps.*sinh(zz./Eratio).*((R-r.^2)+(4*R-r2.^2)-x0*cos(theta).*(R-r.^2)./r)./sinh(1/Eratio);
figure
pcolor(yy,zz,squeeze(vtheta4))
colormap(darkb2r(-1,1))
shading 'flat'
colorbar
axis equal
axis tight

E=0.25;
x0=-0.5;
eps=0.01;

Eratio=sqrt(E);
s=sin(0.5/Eratio);
c=cos(0.5/Eratio);
s2=sinh(0.5/Eratio);
c2=cosh(0.5/Eratio);
A=-0.5*c*s2/((s^2)*(c2^2)+(c^2)*(s2^2));
B=0.5*c2*s/((s^2)*(c2^2)+(c^2)*(s2^2));
D=A*(s*s2-c*c2)+B*(c*c2+s*s2);
r2=sqrt((yy.^2)+(xx-x0).^2);
zeta=(zz-0.5)/Eratio;
vtheta=r.*((R-r).^2).*(0.5+B.*sin(zeta).*cosh(zeta)-A.*cos(zeta).*sinh(zeta));
vtheta5=5.*vtheta-2*eps.*sinh(zz./Eratio).*((R-r.^2)+(4*R-r2.^2)-x0*cos(theta).*(R-r.^2)./r)./sinh(1/Eratio);
figure
pcolor(yy,zz,squeeze(vtheta5))
colormap(darkb2r(-1,1))
shading 'flat'
colorbar
axis equal
axis tight
title('E 0.25, eps 0.01 current')

E=0.25;
x0=-0.5;
eps=0.025;

Eratio=sqrt(E);
s=sin(0.5/Eratio);
c=cos(0.5/Eratio);
s2=sinh(0.5/Eratio);
c2=cosh(0.5/Eratio);
A=-0.5*c*s2/((s^2)*(c2^2)+(c^2)*(s2^2));
B=0.5*c2*s/((s^2)*(c2^2)+(c^2)*(s2^2));
D=A*(s*s2-c*c2)+B*(c*c2+s*s2);
r2=sqrt((yy.^2)+(xx-x0).^2);
zeta=(zz-0.5)/Eratio;
vtheta=r.*((R-r).^2).*(0.5+B.*sin(zeta).*cosh(zeta)-A.*cos(zeta).*sinh(zeta));
vtheta16=5.*vtheta-2*eps.*sinh(zz./Eratio).*((R-r.^2)+(4*R-r2.^2)-x0*cos(theta).*(R-r.^2)./r)./sinh(1/Eratio);

E=0.25;
x0=-0.5;
eps=0.005;

Eratio=sqrt(E);
s=sin(0.5/Eratio);
c=cos(0.5/Eratio);
s2=sinh(0.5/Eratio);
c2=cosh(0.5/Eratio);
A=-0.5*c*s2/((s^2)*(c2^2)+(c^2)*(s2^2));
B=0.5*c2*s/((s^2)*(c2^2)+(c^2)*(s2^2));
D=A*(s*s2-c*c2)+B*(c*c2+s*s2);
r2=sqrt((yy.^2)+(xx-x0).^2);
zeta=(zz-0.5)/Eratio;
vtheta=r.*((R-r).^2).*(0.5+B.*sin(zeta).*cosh(zeta)-A.*cos(zeta).*sinh(zeta));
vtheta17=5.*vtheta-2*eps.*sinh(zz./Eratio).*((R-r.^2)+(4*R-r2.^2)-x0*cos(theta).*(R-r.^2)./r)./sinh(1/Eratio);

E=0.25;
x0=-0.5;
eps=0.0025;

Eratio=sqrt(E);
s=sin(0.5/Eratio);
c=cos(0.5/Eratio);
s2=sinh(0.5/Eratio);
c2=cosh(0.5/Eratio);
A=-0.5*c*s2/((s^2)*(c2^2)+(c^2)*(s2^2));
B=0.5*c2*s/((s^2)*(c2^2)+(c^2)*(s2^2));
D=A*(s*s2-c*c2)+B*(c*c2+s*s2);
r2=sqrt((yy.^2)+(xx-x0).^2);
zeta=(zz-0.5)/Eratio;
vtheta=r.*((R-r).^2).*(0.5+B.*sin(zeta).*cosh(zeta)-A.*cos(zeta).*sinh(zeta));
vtheta18=5.*vtheta-2*eps.*sinh(zz./Eratio).*((R-r.^2)+(4*R-r2.^2)-x0*cos(theta).*(R-r.^2)./r)./sinh(1/Eratio);



E=0.125;
x0=-0.9;
eps=0.05;

Eratio=sqrt(E);
s=sin(0.5/Eratio);
c=cos(0.5/Eratio);
s2=sinh(0.5/Eratio);
c2=cosh(0.5/Eratio);
A=-0.5*c*s2/((s^2)*(c2^2)+(c^2)*(s2^2));
B=0.5*c2*s/((s^2)*(c2^2)+(c^2)*(s2^2));
D=A*(s*s2-c*c2)+B*(c*c2+s*s2);
r2=sqrt((yy.^2)+(xx-x0).^2);
zeta=(zz-0.5)/Eratio;
vtheta=r.*((R-r).^2).*(0.5+B.*sin(zeta).*cosh(zeta)-A.*cos(zeta).*sinh(zeta));
vtheta6=5.*vtheta-2*eps.*sinh(zz./Eratio).*((R-r.^2)+(4*R-r2.^2)-x0*cos(theta).*(R-r.^2)./r)./sinh(1/Eratio);
figure
pcolor(yy,zz,squeeze(vtheta6))
colormap(darkb2r(-1,1))
shading 'flat'
colorbar
axis equal
axis tight
title('E 0.125, eps 0.05 current')

E=0.125;
x0=-0.9;
eps=0.01;

Eratio=sqrt(E);
s=sin(0.5/Eratio);
c=cos(0.5/Eratio);
s2=sinh(0.5/Eratio);
c2=cosh(0.5/Eratio);
A=-0.5*c*s2/((s^2)*(c2^2)+(c^2)*(s2^2));
B=0.5*c2*s/((s^2)*(c2^2)+(c^2)*(s2^2));
D=A*(s*s2-c*c2)+B*(c*c2+s*s2);
r2=sqrt((yy.^2)+(xx-x0).^2);
zeta=(zz-0.5)/Eratio;
vtheta=r.*((R-r).^2).*(0.5+B.*sin(zeta).*cosh(zeta)-A.*cos(zeta).*sinh(zeta));
vtheta7=5.*vtheta-2*eps.*sinh(zz./Eratio).*((R-r.^2)+(4*R-r2.^2)-x0*cos(theta).*(R-r.^2)./r)./sinh(1/Eratio);
figure
pcolor(yy,zz,squeeze(vtheta7))
colormap(darkb2r(-1,1))
shading 'flat'
colorbar
axis equal
axis tight
title('E 0.125, eps 0.01')

E=0.125;
x0=-0.9;
eps=0.1;

Eratio=sqrt(E);
s=sin(0.5/Eratio);
c=cos(0.5/Eratio);
s2=sinh(0.5/Eratio);
c2=cosh(0.5/Eratio);
A=-0.5*c*s2/((s^2)*(c2^2)+(c^2)*(s2^2));
B=0.5*c2*s/((s^2)*(c2^2)+(c^2)*(s2^2));
D=A*(s*s2-c*c2)+B*(c*c2+s*s2);
r2=sqrt((yy.^2)+(xx-x0).^2);
zeta=(zz-0.5)/Eratio;
vtheta=r.*((R-r).^2).*(0.5+B.*sin(zeta).*cosh(zeta)-A.*cos(zeta).*sinh(zeta));
vtheta19=5.*vtheta-2*eps.*sinh(zz./Eratio).*((R-r.^2)+(4*R-r2.^2)-x0*cos(theta).*(R-r.^2)./r)./sinh(1/Eratio);

E=0.125;
x0=-0.9;
eps=0.025;

Eratio=sqrt(E);
s=sin(0.5/Eratio);
c=cos(0.5/Eratio);
s2=sinh(0.5/Eratio);
c2=cosh(0.5/Eratio);
A=-0.5*c*s2/((s^2)*(c2^2)+(c^2)*(s2^2));
B=0.5*c2*s/((s^2)*(c2^2)+(c^2)*(s2^2));
D=A*(s*s2-c*c2)+B*(c*c2+s*s2);
r2=sqrt((yy.^2)+(xx-x0).^2);
zeta=(zz-0.5)/Eratio;
vtheta=r.*((R-r).^2).*(0.5+B.*sin(zeta).*cosh(zeta)-A.*cos(zeta).*sinh(zeta));
vtheta20=5.*vtheta-2*eps.*sinh(zz./Eratio).*((R-r.^2)+(4*R-r2.^2)-x0*cos(theta).*(R-r.^2)./r)./sinh(1/Eratio);

E=0.125;
x0=-0.9;
eps=0.005;

Eratio=sqrt(E);
s=sin(0.5/Eratio);
c=cos(0.5/Eratio);
s2=sinh(0.5/Eratio);
c2=cosh(0.5/Eratio);
A=-0.5*c*s2/((s^2)*(c2^2)+(c^2)*(s2^2));
B=0.5*c2*s/((s^2)*(c2^2)+(c^2)*(s2^2));
D=A*(s*s2-c*c2)+B*(c*c2+s*s2);
r2=sqrt((yy.^2)+(xx-x0).^2);
zeta=(zz-0.5)/Eratio;
vtheta=r.*((R-r).^2).*(0.5+B.*sin(zeta).*cosh(zeta)-A.*cos(zeta).*sinh(zeta));
vtheta21=5.*vtheta-2*eps.*sinh(zz./Eratio).*((R-r.^2)+(4*R-r2.^2)-x0*cos(theta).*(R-r.^2)./r)./sinh(1/Eratio);

E=0.125;
x0=-0.9;
eps=0.0025;

Eratio=sqrt(E);
s=sin(0.5/Eratio);
c=cos(0.5/Eratio);
s2=sinh(0.5/Eratio);
c2=cosh(0.5/Eratio);
A=-0.5*c*s2/((s^2)*(c2^2)+(c^2)*(s2^2));
B=0.5*c2*s/((s^2)*(c2^2)+(c^2)*(s2^2));
D=A*(s*s2-c*c2)+B*(c*c2+s*s2);
r2=sqrt((yy.^2)+(xx-x0).^2);
zeta=(zz-0.5)/Eratio;
vtheta=r.*((R-r).^2).*(0.5+B.*sin(zeta).*cosh(zeta)-A.*cos(zeta).*sinh(zeta));
vtheta22=5.*vtheta-2*eps.*sinh(zz./Eratio).*((R-r.^2)+(4*R-r2.^2)-x0*cos(theta).*(R-r.^2)./r)./sinh(1/Eratio);


E=0.02;
x0=-0.5;
eps=0.1;

Eratio=sqrt(E);
s=sin(0.5/Eratio);
c=cos(0.5/Eratio);
s2=sinh(0.5/Eratio);
c2=cosh(0.5/Eratio);
A=-0.5*c*s2/((s^2)*(c2^2)+(c^2)*(s2^2));
B=0.5*c2*s/((s^2)*(c2^2)+(c^2)*(s2^2));
D=A*(s*s2-c*c2)+B*(c*c2+s*s2);
r2=sqrt((yy.^2)+(xx-x0).^2);
zeta=(zz-0.5)/Eratio;
vtheta=r.*((R-r).^2).*(0.5+B.*sin(zeta).*cosh(zeta)-A.*cos(zeta).*sinh(zeta));
vtheta8=5.*vtheta-2*eps.*sinh(zz./Eratio).*((R-r.^2)+(4*R-r2.^2)-x0*cos(theta).*(R-r.^2)./r)./sinh(1/Eratio);
figure
pcolor(yy,zz,squeeze(vtheta8))
colormap(darkb2r(-1,1))
shading 'flat'
colorbar
axis equal
axis tight
title('E 0.02, eps 0.1 current')

E=0.02;
x0=-0.5;
eps=0.01;

Eratio=sqrt(E);
s=sin(0.5/Eratio);
c=cos(0.5/Eratio);
s2=sinh(0.5/Eratio);
c2=cosh(0.5/Eratio);
A=-0.5*c*s2/((s^2)*(c2^2)+(c^2)*(s2^2));
B=0.5*c2*s/((s^2)*(c2^2)+(c^2)*(s2^2));
D=A*(s*s2-c*c2)+B*(c*c2+s*s2);
r2=sqrt((yy.^2)+(xx-x0).^2);
zeta=(zz-0.5)/Eratio;
vtheta=r.*((R-r).^2).*(0.5+B.*sin(zeta).*cosh(zeta)-A.*cos(zeta).*sinh(zeta));
vtheta9=5.*vtheta-2*eps.*sinh(zz./Eratio).*((R-r.^2)+(4*R-r2.^2)-x0*cos(theta).*(R-r.^2)./r)./sinh(1/Eratio);
figure
pcolor(yy,zz,squeeze(vtheta9))
colormap(darkb2r(-1,1))
shading 'flat'
colorbar
axis equal
axis tight
title('E 0.02, eps 0.01')

E=0.02;
x0=-0.5;
eps=0.05;

Eratio=sqrt(E);
s=sin(0.5/Eratio);
c=cos(0.5/Eratio);
s2=sinh(0.5/Eratio);
c2=cosh(0.5/Eratio);
A=-0.5*c*s2/((s^2)*(c2^2)+(c^2)*(s2^2));
B=0.5*c2*s/((s^2)*(c2^2)+(c^2)*(s2^2));
D=A*(s*s2-c*c2)+B*(c*c2+s*s2);
r2=sqrt((yy.^2)+(xx-x0).^2);
zeta=(zz-0.5)/Eratio;
vtheta=r.*((R-r).^2).*(0.5+B.*sin(zeta).*cosh(zeta)-A.*cos(zeta).*sinh(zeta));
vtheta23=5.*vtheta-2*eps.*sinh(zz./Eratio).*((R-r.^2)+(4*R-r2.^2)-x0*cos(theta).*(R-r.^2)./r)./sinh(1/Eratio);

E=0.02;
x0=-0.5;
eps=0.025;

Eratio=sqrt(E);
s=sin(0.5/Eratio);
c=cos(0.5/Eratio);
s2=sinh(0.5/Eratio);
c2=cosh(0.5/Eratio);
A=-0.5*c*s2/((s^2)*(c2^2)+(c^2)*(s2^2));
B=0.5*c2*s/((s^2)*(c2^2)+(c^2)*(s2^2));
D=A*(s*s2-c*c2)+B*(c*c2+s*s2);
r2=sqrt((yy.^2)+(xx-x0).^2);
zeta=(zz-0.5)/Eratio;
vtheta=r.*((R-r).^2).*(0.5+B.*sin(zeta).*cosh(zeta)-A.*cos(zeta).*sinh(zeta));
vtheta24=5.*vtheta-2*eps.*sinh(zz./Eratio).*((R-r.^2)+(4*R-r2.^2)-x0*cos(theta).*(R-r.^2)./r)./sinh(1/Eratio);

E=0.02;
x0=-0.5;
eps=0.005;

Eratio=sqrt(E);
s=sin(0.5/Eratio);
c=cos(0.5/Eratio);
s2=sinh(0.5/Eratio);
c2=cosh(0.5/Eratio);
A=-0.5*c*s2/((s^2)*(c2^2)+(c^2)*(s2^2));
B=0.5*c2*s/((s^2)*(c2^2)+(c^2)*(s2^2));
D=A*(s*s2-c*c2)+B*(c*c2+s*s2);
r2=sqrt((yy.^2)+(xx-x0).^2);
zeta=(zz-0.5)/Eratio;
vtheta=r.*((R-r).^2).*(0.5+B.*sin(zeta).*cosh(zeta)-A.*cos(zeta).*sinh(zeta));
vtheta25=5.*vtheta-2*eps.*sinh(zz./Eratio).*((R-r.^2)+(4*R-r2.^2)-x0*cos(theta).*(R-r.^2)./r)./sinh(1/Eratio);

E=0.02;
x0=-0.5;
eps=0.0025;

Eratio=sqrt(E);
s=sin(0.5/Eratio);
c=cos(0.5/Eratio);
s2=sinh(0.5/Eratio);
c2=cosh(0.5/Eratio);
A=-0.5*c*s2/((s^2)*(c2^2)+(c^2)*(s2^2));
B=0.5*c2*s/((s^2)*(c2^2)+(c^2)*(s2^2));
D=A*(s*s2-c*c2)+B*(c*c2+s*s2);
r2=sqrt((yy.^2)+(xx-x0).^2);
zeta=(zz-0.5)/Eratio;
vtheta=r.*((R-r).^2).*(0.5+B.*sin(zeta).*cosh(zeta)-A.*cos(zeta).*sinh(zeta));
vtheta26=5.*vtheta-2*eps.*sinh(zz./Eratio).*((R-r.^2)+(4*R-r2.^2)-x0*cos(theta).*(R-r.^2)./r)./sinh(1/Eratio);

E=0.0005;
x0=-0.9;
eps=0.01;

Eratio=sqrt(E);
s=sin(0.5/Eratio);
c=cos(0.5/Eratio);
s2=sinh(0.5/Eratio);
c2=cosh(0.5/Eratio);
A=-0.5*c*s2/((s^2)*(c2^2)+(c^2)*(s2^2));
B=0.5*c2*s/((s^2)*(c2^2)+(c^2)*(s2^2));
D=A*(s*s2-c*c2)+B*(c*c2+s*s2);
r2=sqrt((yy.^2)+(xx-x0).^2);
zeta=(zz-0.5)/Eratio;
vtheta=r.*((R-r).^2).*(0.5+B.*sin(zeta).*cosh(zeta)-A.*cos(zeta).*sinh(zeta));
vtheta10=5.*vtheta-2*eps.*sinh(zz./Eratio).*((R-r.^2)+(4*R-r2.^2)-x0*cos(theta).*(R-r.^2)./r)./sinh(1/Eratio);
figure
pcolor(yy,zz,squeeze(vtheta10))
colormap(darkb2r(-1,1))
shading 'flat'
colorbar
axis equal
axis tight
title('E 0.0005, eps 0.01 current')

E=0.0005;
x0=-0.9;
eps=0.05;

Eratio=sqrt(E);
s=sin(0.5/Eratio);
c=cos(0.5/Eratio);
s2=sinh(0.5/Eratio);
c2=cosh(0.5/Eratio);
A=-0.5*c*s2/((s^2)*(c2^2)+(c^2)*(s2^2));
B=0.5*c2*s/((s^2)*(c2^2)+(c^2)*(s2^2));
D=A*(s*s2-c*c2)+B*(c*c2+s*s2);
r2=sqrt((yy.^2)+(xx-x0).^2);
zeta=(zz-0.5)/Eratio;
vtheta=r.*((R-r).^2).*(0.5+B.*sin(zeta).*cosh(zeta)-A.*cos(zeta).*sinh(zeta));
vtheta11=5.*vtheta-2*eps.*sinh(zz./Eratio).*((R-r.^2)+(4*R-r2.^2)-x0*cos(theta).*(R-r.^2)./r)./sinh(1/Eratio);
figure
pcolor(yy,zz,squeeze(vtheta11))
colormap(darkb2r(-1,1))
shading 'flat'
colorbar
axis equal
axis tight
title('E 0.0005, eps 0.05')

E=0.0005;
x0=-0.9;
eps=0.1;

Eratio=sqrt(E);
s=sin(0.5/Eratio);
c=cos(0.5/Eratio);
s2=sinh(0.5/Eratio);
c2=cosh(0.5/Eratio);
A=-0.5*c*s2/((s^2)*(c2^2)+(c^2)*(s2^2));
B=0.5*c2*s/((s^2)*(c2^2)+(c^2)*(s2^2));
D=A*(s*s2-c*c2)+B*(c*c2+s*s2);
r2=sqrt((yy.^2)+(xx-x0).^2);
zeta=(zz-0.5)/Eratio;
vtheta=r.*((R-r).^2).*(0.5+B.*sin(zeta).*cosh(zeta)-A.*cos(zeta).*sinh(zeta));
vtheta12=5.*vtheta-2*eps.*sinh(zz./Eratio).*((R-r.^2)+(4*R-r2.^2)-x0*cos(theta).*(R-r.^2)./r)./sinh(1/Eratio);

E=0.0005;
x0=-0.9;
eps=0.025;

Eratio=sqrt(E);
s=sin(0.5/Eratio);
c=cos(0.5/Eratio);
s2=sinh(0.5/Eratio);
c2=cosh(0.5/Eratio);
A=-0.5*c*s2/((s^2)*(c2^2)+(c^2)*(s2^2));
B=0.5*c2*s/((s^2)*(c2^2)+(c^2)*(s2^2));
D=A*(s*s2-c*c2)+B*(c*c2+s*s2);
r2=sqrt((yy.^2)+(xx-x0).^2);
zeta=(zz-0.5)/Eratio;
vtheta=r.*((R-r).^2).*(0.5+B.*sin(zeta).*cosh(zeta)-A.*cos(zeta).*sinh(zeta));
vtheta13=5.*vtheta-2*eps.*sinh(zz./Eratio).*((R-r.^2)+(4*R-r2.^2)-x0*cos(theta).*(R-r.^2)./r)./sinh(1/Eratio);

E=0.0005;
x0=-0.9;
eps=0.005;

Eratio=sqrt(E);
s=sin(0.5/Eratio);
c=cos(0.5/Eratio);
s2=sinh(0.5/Eratio);
c2=cosh(0.5/Eratio);
A=-0.5*c*s2/((s^2)*(c2^2)+(c^2)*(s2^2));
B=0.5*c2*s/((s^2)*(c2^2)+(c^2)*(s2^2));
D=A*(s*s2-c*c2)+B*(c*c2+s*s2);
r2=sqrt((yy.^2)+(xx-x0).^2);
zeta=(zz-0.5)/Eratio;
vtheta=r.*((R-r).^2).*(0.5+B.*sin(zeta).*cosh(zeta)-A.*cos(zeta).*sinh(zeta));
vtheta14=5.*vtheta-2*eps.*sinh(zz./Eratio).*((R-r.^2)+(4*R-r2.^2)-x0*cos(theta).*(R-r.^2)./r)./sinh(1/Eratio);

E=0.0005;
x0=-0.9;
eps=0.0025;

Eratio=sqrt(E);
s=sin(0.5/Eratio);
c=cos(0.5/Eratio);
s2=sinh(0.5/Eratio);
c2=cosh(0.5/Eratio);
A=-0.5*c*s2/((s^2)*(c2^2)+(c^2)*(s2^2));
B=0.5*c2*s/((s^2)*(c2^2)+(c^2)*(s2^2));
D=A*(s*s2-c*c2)+B*(c*c2+s*s2);
r2=sqrt((yy.^2)+(xx-x0).^2);
zeta=(zz-0.5)/Eratio;
vtheta=r.*((R-r).^2).*(0.5+B.*sin(zeta).*cosh(zeta)-A.*cos(zeta).*sinh(zeta));
vtheta15=5.*vtheta-2*eps.*sinh(zz./Eratio).*((R-r.^2)+(4*R-r2.^2)-x0*cos(theta).*(R-r.^2)./r)./sinh(1/Eratio);


figure
subplot(6,4,1)
pcolor(squeeze(yy),squeeze(zz),squeeze(vtheta2))
colormap(darkb2r(-1,1))
shading 'flat'
axis equal
axis tight
hold on
contour(squeeze(yy),squeeze(zz),squeeze(vtheta2),[0,0],'k')
set(gca,'XTickLabel',[]);
ylabel('epsilon 0.1')
subplot(6,4,2)
pcolor(squeeze(yy),squeeze(zz),squeeze(vtheta19))
colormap(darkb2r(-1,1))
shading 'flat'
axis equal
axis tight
hold on
contour(squeeze(yy),squeeze(zz),squeeze(vtheta19),[0,0],'k')
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
subplot(6,4,3)
pcolor(squeeze(yy),squeeze(zz),squeeze(vtheta8))
colormap(darkb2r(-1,1))
shading 'flat'
axis equal
axis tight
hold on
contour(squeeze(yy),squeeze(zz),squeeze(vtheta8),[0,0],'k')
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
subplot(6,4,4)
pcolor(squeeze(yy),squeeze(zz),squeeze(vtheta12))
colormap(darkb2r(-1,1))
shading 'flat'
axis equal
axis tight
hold on
contour(squeeze(yy),squeeze(zz),squeeze(vtheta12),[0,0],'k')
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
subplot(6,4,5)
pcolor(squeeze(yy),squeeze(zz),squeeze(vtheta4))
colormap(darkb2r(-1,1))
shading 'flat'
axis equal
axis tight
hold on
contour(squeeze(yy),squeeze(zz),squeeze(vtheta4),[0,0],'k')
set(gca,'XTickLabel',[]);
ylabel('epsilon 0.05')
subplot(6,4,6)
pcolor(squeeze(yy),squeeze(zz),squeeze(vtheta6))
colormap(darkb2r(-1,1))
shading 'flat'
axis equal
axis tight
hold on
contour(squeeze(yy),squeeze(zz),squeeze(vtheta6),[0,0],'k')
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
subplot(6,4,7)
pcolor(squeeze(yy),squeeze(zz),squeeze(vtheta23))
colormap(darkb2r(-1,1))
shading 'flat'
axis equal
axis tight
hold on
contour(squeeze(yy),squeeze(zz),squeeze(vtheta23),[0,0],'k')
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
subplot(6,4,8)
pcolor(squeeze(yy),squeeze(zz),squeeze(vtheta11))
colormap(darkb2r(-1,1))
shading 'flat'
axis equal
axis tight
hold on
contour(squeeze(yy),squeeze(zz),squeeze(vtheta11),[0,0],'k')
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
subplot(6,4,9)
pcolor(squeeze(yy),squeeze(zz),squeeze(vtheta16))
colormap(darkb2r(-1,1))
shading 'flat'
axis equal
axis tight
hold on
contour(squeeze(yy),squeeze(zz),squeeze(vtheta16),[0,0],'k')
set(gca,'XTickLabel',[]);
ylabel('epsilon 0.025')
subplot(6,4,10)
pcolor(squeeze(yy),squeeze(zz),squeeze(vtheta20))
colormap(darkb2r(-1,1))
shading 'flat'
axis equal
axis tight
hold on
contour(squeeze(yy),squeeze(zz),squeeze(vtheta20),[0,0],'k')
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
subplot(6,4,11)
pcolor(squeeze(yy),squeeze(zz),squeeze(vtheta24))
colormap(darkb2r(-1,1))
shading 'flat'
axis equal
axis tight
hold on
contour(squeeze(yy),squeeze(zz),squeeze(vtheta24),[0,0],'k')
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
subplot(6,4,12)
pcolor(squeeze(yy),squeeze(zz),squeeze(vtheta13))
colormap(darkb2r(-1,1))
shading 'flat'
axis equal
axis tight
hold on
contour(squeeze(yy),squeeze(zz),squeeze(vtheta13),[0,0],'k')
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
subplot(6,4,13)
pcolor(squeeze(yy),squeeze(zz),squeeze(vtheta5))
colormap(darkb2r(-1,1))
shading 'flat'
axis equal
axis tight
hold on
contour(squeeze(yy),squeeze(zz),squeeze(vtheta5),[0,0],'k')
set(gca,'XTickLabel',[]);
ylabel('epsilon 0.01')
subplot(6,4,14)
pcolor(squeeze(yy),squeeze(zz),squeeze(vtheta7))
colormap(darkb2r(-1,1))
shading 'flat'
axis equal
axis tight
hold on
contour(squeeze(yy),squeeze(zz),squeeze(vtheta7),[0,0],'k')
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
subplot(6,4,15)
pcolor(squeeze(yy),squeeze(zz),squeeze(vtheta9))
colormap(darkb2r(-1,1))
shading 'flat'
axis equal
axis tight
hold on
contour(squeeze(yy),squeeze(zz),squeeze(vtheta9),[0,0],'k')
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
subplot(6,4,16)
pcolor(squeeze(yy),squeeze(zz),squeeze(vtheta10))
colormap(darkb2r(-1,1))
shading 'flat'
axis equal
axis tight
hold on
contour(squeeze(yy),squeeze(zz),squeeze(vtheta10),[0,0],'k')
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
subplot(6,4,17)
pcolor(squeeze(yy),squeeze(zz),squeeze(vtheta17))
colormap(darkb2r(-1,1))
shading 'flat'
axis equal
axis tight
hold on
contour(squeeze(yy),squeeze(zz),squeeze(vtheta17),[0,0],'k')
set(gca,'XTickLabel',[]);
ylabel('epsilon 0.005')
subplot(6,4,18)
pcolor(squeeze(yy),squeeze(zz),squeeze(vtheta21))
colormap(darkb2r(-1,1))
shading 'flat'
axis equal
axis tight
hold on
contour(squeeze(yy),squeeze(zz),squeeze(vtheta21),[0,0],'k')
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
subplot(6,4,19)
pcolor(squeeze(yy),squeeze(zz),squeeze(vtheta25))
colormap(darkb2r(-1,1))
shading 'flat'
axis equal
axis tight
hold on
contour(squeeze(yy),squeeze(zz),squeeze(vtheta25),[0,0],'k')
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
subplot(6,4,20)
pcolor(squeeze(yy),squeeze(zz),squeeze(vtheta14))
colormap(darkb2r(-1,1))
shading 'flat'
axis equal
axis tight
hold on
contour(squeeze(yy),squeeze(zz),squeeze(vtheta14),[0,0],'k')
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
subplot(6,4,21)
pcolor(squeeze(yy),squeeze(zz),squeeze(vtheta18))
colormap(darkb2r(-1,1))
shading 'flat'
axis equal
axis tight
hold on
contour(squeeze(yy),squeeze(zz),squeeze(vtheta18),[0,0],'k')
ylabel('epsilon 0.0025')
xlabel('E 0.25')
subplot(6,4,22)
pcolor(squeeze(yy),squeeze(zz),squeeze(vtheta22))
colormap(darkb2r(-1,1))
shading 'flat'
axis equal
axis tight
hold on
contour(squeeze(yy),squeeze(zz),squeeze(vtheta22),[0,0],'k')
set(gca,'YTickLabel',[]);
xlabel('E 0.125')
subplot(6,4,23)
pcolor(squeeze(yy),squeeze(zz),squeeze(vtheta26))
colormap(darkb2r(-1,1))
shading 'flat'
axis equal
axis tight
hold on
contour(squeeze(yy),squeeze(zz),squeeze(vtheta26),[0,0],'k')
set(gca,'YTickLabel',[]);
xlabel('E 0.02')
subplot(6,4,24)
pcolor(squeeze(yy),squeeze(zz),squeeze(vtheta15))
colormap(darkb2r(-1,1))
shading 'flat'
axis equal
axis tight
hold on
contour(squeeze(yy),squeeze(zz),squeeze(vtheta15),[0,0],'k')
set(gca,'YTickLabel',[]);
xlabel('E 0.0005')

%%
%using eps 0.01, make poincare sections
load('E0.25P4e0.01x0-0.5.mat')
figure
hold all;
for itr=1:size(ztr,2);%[15 59 70 85];%:size(xtr,2)
             jj=find(xtr(1:end-1,itr).*xtr(2:end,itr)<0); % look for 
%points along the trajectory where x changes sign
             plot(ytr(jj,itr),ztr(jj,itr),'k.','MarkerSize',1);
end
xlabel('Y','fontsize',14)
 ylabel('Z','fontsize',14)
 title(strcat('E=',num2str(E)),'fontsize',16)
 set(gca,'xtick',[-1 0 1])
 set(gca,'ytick',[0 0.5 1])
 axis equal
axis tight

figure
load('E0.125P4e0.01x0-0.9.mat')
hold all;
for itr=1:size(ztr,2);%[15 59 70 85];%:size(xtr,2)
             jj=find(xtr(1:end-1,itr).*xtr(2:end,itr)<0); % look for 
%points along the trajectory where x changes sign
             plot(ytr(jj,itr),ztr(jj,itr),'k.','MarkerSize',1);

end
 title(strcat('E=',num2str(E)),'fontsize',16)
 set(gca,'xtick',[])
 set(gca,'ytick',[])
 axis equal
axis tight

figure
load('E0.02P4e0.01x0-0.5.mat')
hold all;
for itr=1:size(ztr,2);%[15 59 70 85];%:size(xtr,2)
             jj=find(xtr(1:end-1,itr).*xtr(2:end,itr)<0); % look for 
%points along the trajectory where x changes sign
             plot(ytr(jj,itr),ztr(jj,itr),'k.','MarkerSize',1);
end
 title(strcat('E=',num2str(E)),'fontsize',16)
 set(gca,'xtick',[])
 set(gca,'ytick',[])
 axis equal
axis tight

figure
load('E0.0005P4e0.01x0-0.9.mat')
hold all;
for itr=1:size(ztr,2);%[15 59 70 85];%:size(xtr,2)
             jj=find(xtr(1:end-1,itr).*xtr(2:end,itr)<0); % look for 
%points along the trajectory where x changes sign
             plot(ytr(jj,itr),ztr(jj,itr),'k.','MarkerSize',1);

end
 title(strcat('E=',num2str(E)),'fontsize',16)
 set(gca,'xtick',[])
 set(gca,'ytick',[])
 axis equal
axis tight
%%
%Plot eps=0.01 poincare on the matching FTLE

%E=0.25
figure
load('ftleE025.mat')
i=4;
[~,Fp,Fn,Fm,~]=ftleCalc(xin3b,yin3b,zin3b,xtr(i+1,:),ytr(i+1,:),ztr(i+1,:),i);
%subplot('position',[0.16 0.1 0.1 0.8])
pcolor(yin3b(2:end-1,1,2:end-1),zin3b(2:end-1,1,2:end-1),squeeze(Fp))
%shading flat
axis equal
axis tight
colormap(spring)
colorbar('SouthOutside')
title('E=0.25','fontsize',16)
xlabel('y','fontsize',14)
ylabel('z','fontsize',14)
set(gca,'xtick',[0.25 0.75])
set(gca,'ytick',[0.25 0.75])
load('E0.25P4e0.01x0-0.5.mat')
%figure
hold all;
for itr=1:size(ztr,2);%[15 59 70 85];%:size(xtr,2)
             jj=find(xtr(1:end-1,itr).*xtr(2:end,itr)<0); % look for 
%points along the trajectory where x changes sign
             plot(ytr(jj,itr),ztr(jj,itr),'k.','MarkerSize',1);
end
xlabel('Y','fontsize',14)
 ylabel('Z','fontsize',14)
 title(strcat('E=',num2str(E)),'fontsize',16)
 set(gca,'xtick',[-1 0 1])
 set(gca,'ytick',[0 0.5 1])
 axis equal
axis tight


%E=0.0005
load('ftleP4E00005.mat')
i=5;
[~,Fp,Fn,Fm,~]=ftleCalc(xin3b,yin3b,zin3b,xtr(i+1,:),ytr(i+1,:),ztr(i+1,:),i);
figure
pcolor(yin3b(2:end-1,1,2:end-1),zin3b(2:end-1,1,2:end-1),squeeze(Fp))
%shading flat
axis equal
axis tight
colormap(spring)
colorbar('SouthOutside')
title('E=0.0005','fontsize',16)
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
%figure
load('E0.0005P4e0.01x0-0.9.mat')
hold all;
for itr=1:size(ztr,2);%[15 59 70 85];%:size(xtr,2)
             jj=find(xtr(1:end-1,itr).*xtr(2:end,itr)<0); % look for 
%points along the trajectory where x changes sign
             plot(ytr(jj,itr),ztr(jj,itr),'k.','MarkerSize',1);

end
 title(strcat('E=',num2str(E)),'fontsize',16)
 set(gca,'xtick',[])
 set(gca,'ytick',[])
 axis equal
axis tight

%%
%Plotting pert4

yin=-1:0.01:1;
xin=-1:0.01:1;
[xx,yy,zz]=meshgrid(xin,yin,1);
R=1;
[theta,r]=cart2pol(xx,yy);

E=0.25;
x0=-0.5;
Eratio=sqrt(E);
r2=sqrt((xx-x0).^2 +yy.^2);

psi1=-sinh(zz/Eratio).*(R.^2 -r.^2).*(4*R.^2 -r2.^2)./sinh(1/Eratio);
psi2=-sinh(zz/Eratio).*(R.^2 -r.^2).*(1*R.^2 -r2.^2)./sinh(1/Eratio);
psi3=-sinh(zz/Eratio).*(R.^2 -r.^2).*(9*R.^2 -r2.^2)./sinh(1/Eratio);

psi1(r>R)=NaN;
psi2(r>R)=NaN;
psi3(r>R)=NaN;

figure
subplot(1,3,1)
contourf(squeeze(xx),squeeze(yy),squeeze(psi2),20)
colormap(darkb2r(-8,1))
axis equal
axis tight
set(gca,'color',[0.5 0.5 0.5])
subplot(1,3,2)
contourf(squeeze(xx),squeeze(yy),squeeze(psi1),20)
colormap(darkb2r(-8,1))
axis equal
axis tight
set(gca,'color',[0.5 0.5 0.5])
subplot(1,3,3)
contourf(squeeze(xx),squeeze(yy),squeeze(psi3),20)
colormap(darkb2r(-8,1))
axis equal
axis tight
colorbar('EastOutside')
set(gca,'color',[0.5 0.5 0.5])

yin=-1:0.01:1;
zin=0:0.01:1;
xin=yin;
[xx,yy,zz]=meshgrid(0,yin,zin);
R=1;
[theta,r]=cart2pol(xx,yy);
[xb,yb,zb]=meshgrid(xin,0,zin);
[thetab,rb]=cart2pol(xb,yb);

E=0.125;
x0=-0.5;
Eratio=sqrt(E);
r2=sqrt((xx-x0).^2 +yy.^2);
r2b=sqrt((xb-x0).^2 +yb.^2);

s=sin(0.5/Eratio);
c=cos(0.5/Eratio);
s2=sinh(0.5/Eratio);
c2=cosh(0.5/Eratio);
A=-0.5*c*s2/((s^2)*(c2^2)+(c^2)*(s2^2));
B=0.5*c2*s/((s^2)*(c2^2)+(c^2)*(s2^2));
D=A*(s*s2-c*c2)+B*(c*c2+s*s2);
zeta=(zz-0.5)/Eratio;
vthetaP=-2*eps.*sinh(zz./Eratio).*((R-r.^2)+(4*R-r2.^2)-x0*cos(theta).*(R-r.^2)./r)./sinh(1/Eratio);
urP=2*eps*x0*sinh(zz./Eratio).*sin(theta).*(R-r.^2)./sinh(1/Eratio);

uP=eps.*(2.*yy.*(R - r.^2) + 2.*yy.*(4*R^2 - r2.^2)).*(sinh(zz./Eratio)./sinh(1./Eratio));
vP=-eps.*(2.*(-x0 + xx).*(R - r.^2) + 2.*xx.*(4*R^2 - r2.^2)).*(sinh(zz./Eratio)./sinh(1./Eratio));
urP2=(uP.*cos(theta)+vP.*sin(theta));

vthetaPb=-2*eps.*sinh(zb./Eratio).*((R-rb.^2)+(4*R-r2b.^2)-x0*cos(thetab).*(R-rb.^2)./rb)./sinh(1/Eratio);
urPb=2*eps*x0*sinh(zb./Eratio).*sin(thetab).*(R-rb.^2)./sinh(1/Eratio);

uPb=eps.*(2.*yb.*(R - rb.^2) + 2.*yb.*(4*R^2 - r2b.^2)).*(sinh(zb./Eratio)./sinh(1./Eratio));
vPb=-eps.*(2.*(-x0 + xb).*(R - rb.^2) + 2.*xb.*(4*R^2 - r2b.^2)).*(sinh(zb./Eratio)./sinh(1./Eratio));

figure
subplot(2,2,1)
pcolor(squeeze(yy),squeeze(zz),squeeze(vP));
shading 'flat'
axis equal
axis tight
title('vy perturbation')
xlabel('y')
ylabel('z')
colorbar

subplot(2,2,2)
pcolor(squeeze(yy),squeeze(zz),squeeze(uP));
shading 'flat'
axis equal
axis tight
colorbar
title('ux perturbation')
xlabel('y')

subplot(2,2,3)
pcolor(squeeze(xb),squeeze(zb),squeeze(vPb));
shading 'flat'
axis equal
axis tight
title('vy perturbation')
xlabel('x')
ylabel('z')
colorbar

subplot(2,2,4)
pcolor(squeeze(xb),squeeze(zb),squeeze(uPb));
shading 'flat'
axis equal
axis tight
colorbar
title('ux perturbation')
xlabel('x')

figure
subplot(1,2,1)
pcolor(squeeze(yy),squeeze(zz),squeeze(vthetaP));
shading 'flat'
axis equal
axis tight
title('vtheta perturbation')
colorbar
subplot(1,2,2)
pcolor(squeeze(yy),squeeze(zz),squeeze(urP));
shading 'flat'
axis equal
axis tight
colorbar
title('ur perturbation')

%%
%Plot E=0.02 stochastic to show manifolds AND alignment at long times with
%psi contours

%load('E0.01sigma0005b.mat')
figure
subplot(1,4,1)
r1=r(1:50,:);
z1=ztr(1:50,:);
scatter(r1(:),z1(:),1,'b')
axis([0.4 0.6 0.2 0.8])
title('t=50','fontsize',14)
ylabel('z','fontsize',12)
xlabel('r','fontsize',12)

subplot(1,4,2)
r1=r(1:75,:);
z1=ztr(1:75,:);
scatter(r1(:),z1(:),1,'b')
axis([0.4 0.6 0.2 0.8])
title('t=75','fontsize',14)

subplot(1,4,3)
r1=r(1:150,:);
z1=ztr(1:150,:);
scatter(r1(:),z1(:),1)
axis([0.4 0.6 0.2 0.8])
title('t=150','fontsize',14)

subplot(1,4,4)
r1=r(:);
z1=ztr(:);

x=0.2:0.005:0.8;
z=0.2:0.001:0.8;
[xx,zz]=meshgrid(x,z);
E=0.01;
Eratio=sqrt(E);%ratio of Ekman layer thickness to dimensional tank depth
s=sin(0.5/Eratio);
c=cos(0.5/Eratio);
s2=sinh(0.5/Eratio);
c2=cosh(0.5/Eratio);
A=-0.5*c*s2/((s^2)*(c2^2)+(c^2)*(s2^2));
B=0.5*c2*s/((s^2)*(c2^2)+(c^2)*(s2^2));
D=A*(s*s2-c*c2)+B*(c*c2+s*s2);
%coordinates
rb=abs(xx);
zeta=(zz-0.5)/Eratio;
%streamfunction; 
Rfun=0.5.*(rb.^2).*((R-rb).^2);
Ffun=(A.*(sin(zeta).*sinh(zeta)-cos(zeta).*cosh(zeta))...
    +B.*(cos(zeta).*cosh(zeta)+sin(zeta).*sinh(zeta))-D);
psi=-Eratio.*Rfun.*Ffun;

contourf(xx,zz,psi,[16.2 15.9 15.6 15.3 15 14.2 13 11.5 10 8 6 4 1 0].*1e-4);
colormap(gray)
caxis([-0.0004 max(max(psi))])
hold all
scatter(r1(:),z1(:),1,'b')
axis([0.4 0.6 0.2 0.8])
title('t=10000','fontsize',14)
