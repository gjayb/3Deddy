%plot K_eff figures (early and late cross-sections, full-cylinder K, R_K,
%and chi, and region mean K

% title('\sigma=26 manifold pdf and mean depth','fontsize',24)
% xlabel('longitude','fontsize',24);ylabel('latitude','fontsize',24)
% set(gca,'fontsize',12)

%% new Sept 2018
load('E002compareLowDiffV2.mat')
cLow1=cAllnorm;
cLow0=cAllnorm0;
kLow1=keffmap;
kLow0=keffmap0;
load('E002compareV2.mat')
%% new Sept 2018, paper figure 9
figure;
h1=subplot(4,3,1);
contourf(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(cLow1(50,:,:,5)),[0 0.1 0.3 0.5 0.7 0.9 1])
caxis([0 1]);ylabel({'Dye, t=39';'z'},'fontsize',14); title('x_0=-0.02 k=10^{-6}','fontsize',14);
set(h1,'XTickLabel',{}); 
xlabel('(a)','fontsize',14)
set(gca,'fontsize',14)
colormap(h1,cbrewer('seq','Blues',10))
freezeColors
h2=subplot(4,3,2);
contourf(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(cAllnorm0(50,:,:,5)),0:0.1:1)
caxis([0 1]);  title('x_0=0 k=10^{-4}','fontsize',14);
xlabel('(b)','fontsize',14)
set(gca,'fontsize',14)
set(h2,'XTickLabel',{}); 
set(h2,'YTickLabel',{}); 
colormap(h2,cbrewer('seq','Blues',10))
freezeColors
h3=subplot(4,3,3);
contourf(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(cAllnorm(50,:,:,5)),0:0.1:1)
caxis([0 1]);  title('x_0=-0.02 k=10^{-4}','fontsize',14);
set(gca,'fontsize',14)
set(h3,'XTickLabel',{}); 
set(h3,'YTickLabel',{});
xlabel('(c)','fontsize',14)
cRow1=cbrewer('seq','Blues',10);
freezeColors
colormap(h3,cRow1); cb1=colorbar; cbfreeze(cb1);
 
h4=subplot(4,3,4); 
pcolor(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(kLow1(50,:,:,5)))%,0:0.02:0.12)
caxis([0 0.15]); shading 'flat'
ylabel({'K_{eff}, t=39';'z'},'fontsize',14); 
set(h4,'XTickLabel',{}); 
set(gca,'fontsize',14)
xlabel('(d)','fontsize',14)
colormap(h4,cbrewer('seq','Reds',15))
freezeColors
h5=subplot(4,3,5); 
pcolor(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(keffmap0(50,:,:,5)))%,0:0.02:0.12)
caxis([0 0.15]); shading 'flat'
set(h5,'XTickLabel',{}); 
set(h5,'YTickLabel',{}); 
xlabel('(e)','fontsize',14)
set(gca,'fontsize',14)
colormap(h5,cbrewer('seq','Reds',15))
freezeColors
h6=subplot(4,3,6); 
pcolor(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(keffmap(50,:,:,5)))%,0:0.02:0.12)
caxis([0 0.15]); shading 'flat'
set(h6,'XTickLabel',{}); 
set(h6,'YTickLabel',{}); 
xlabel('(f)','fontsize',14)
set(gca,'fontsize',14)
cRow2=cbrewer('seq','Reds',15);
freezeColors
colormap(cRow2); cb2=colorbar; cbfreeze(cb2);

h7=subplot(4,3,7);
contourf(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(cLow1(50,:,:,31)),[0 0.2:0.05:0.65])
caxis([0.3 0.7]);ylabel({'Dye, t=299';'z'},'fontsize',14); 
set(h7,'XTickLabel',{}); 
set(gca,'fontsize',14)
xlabel('(g)','fontsize',14)
colormap(h7,cbrewer('seq','Blues',40))
freezeColors
h8=subplot(4,3,8);
contourf(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(cAllnorm0(50,:,:,31)),0.455:0.005:0.6)
caxis([0.3 0.7]);  
set(gca,'fontsize',14)
set(h8,'XTickLabel',{}); 
set(h8,'YTickLabel',{}); 
xlabel('(h)','fontsize',14)
colormap(h8,cbrewer('seq','Blues',40))
freezeColors
h9=subplot(4,3,9);
contourf(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(cAllnorm(50,:,:,31)),0.455:0.005:0.6)
caxis([0.3 0.7]); 
set(gca,'fontsize',14)
set(h9,'XTickLabel',{}); 
set(h9,'YTickLabel',{});
xlabel('(i)','fontsize',14)
colormap(h9,cbrewer('seq','Blues',40)); cb3=colorbar; cbfreeze(cb3); freezeColors;

h10=subplot(4,3,10); 
pcolor(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(kLow1(50,:,:,31)))%,0:0.02:0.12)
caxis([0 0.04]); shading 'flat'
ylabel({'K_{eff}, t=299';'z'},'fontsize',14); 
xlabel({'x';'(j)'},'fontsize',14); 
set(gca,'fontsize',14)
colormap(h10,cbrewer('seq','Reds',15))
freezeColors
h11=subplot(4,3,11); 
pcolor(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(keffmap0(50,:,:,31)))%,0:0.02:0.12)
caxis([0 0.04]); shading 'flat'
set(h11,'YTickLabel',{}); 
xlabel({'x';'(k)'},'fontsize',14); 
set(gca,'fontsize',14)
colormap(h11,cbrewer('seq','Reds',15))
freezeColors
h12=subplot(4,3,12); 
pcolor(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(keffmap(50,:,:,31)))%,0:0.02:0.12)
caxis([0 0.04]); shading 'flat'
set(h12,'YTickLabel',{}); 
xlabel({'x';'(l)'},'fontsize',14); 
set(gca,'fontsize',14)
cRow4=cbrewer('seq','Reds',15);
freezeColors
colormap(cRow4); cb4=colorbar; cbfreeze(cb4); 
%% loading for paper figure 10
load('E0125compare3v2.mat')


%% make paper figure 10
figure; 
h1=subplot(4,3,1);
contourf(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(cAllnorm0(50,:,:,5)),0:0.1:1)
caxis([0 1]); 
ylabel({'Dye, t=39';'z'},'fontsize',14); title('x_0=0','fontsize',14); 
set(h1,'XTickLabel',{});
xlabel('(a)','fontsize',14)
set(gca,'fontsize',12)
colormap(cbrewer('seq','Blues',10)); freezeColors;
h2=subplot(4,3,2);
contourf(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(cAllnorm(50,:,:,5)),0:0.1:1)
caxis([0 1]); 
title('x_0=-0.02','fontsize',14); 
set(h2,'XTickLabel',{});
set(h2,'YTickLabel',{});
xlabel('(b)','fontsize',14)
set(gca,'fontsize',12)
colormap(cbrewer('seq','Blues',10)); freezeColors;
h3=subplot(4,3,3);
contourf(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(cAllnorm2(50,:,:,5)),0:0.1:1)
caxis([0 1]); 
title('x_0=-0.02','fontsize',14); 
set(h3,'XTickLabel',{});
xlabel('(c)','fontsize',14)
set(h3,'YTickLabel',{});
set(gca,'fontsize',12)
colormap(cbrewer('seq','Blues',10)); cb1=colorbar; freezeColors; cbfreeze(cb1);

h4=subplot(4,3,4);
pcolor(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(keffmap0(50,:,:,5)))%,0.1:0.02:0.57)
shading 'flat'
caxis([0 0.015]); 
set(gca,'fontsize',12)
ylabel({'K_{eff}, t=39';'z'},'fontsize',14);
set(h4,'XTickLabel',{});
xlabel('(d)','fontsize',14)
colormap(cbrewer('seq','Reds',20)); freezeColors;
h5=subplot(4,3,5);
pcolor(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(keffmap(50,:,:,5)))%,0.1:0.02:0.57)
shading 'flat'
caxis([0 0.015]); 
set(gca,'fontsize',12)
set(h5,'XTickLabel',{});
xlabel('(e)','fontsize',14)
set(h5,'YTickLabel',{});
colormap(cbrewer('seq','Reds',20)); freezeColors;
h6=subplot(4,3,6);
pcolor(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(keffmap2(50,:,:,5)))%,0.1:0.02:0.57)
shading 'flat'
caxis([0 0.015]); 
set(gca,'fontsize',12)
set(h6,'XTickLabel',{});
xlabel('(f)','fontsize',14)
set(h6,'YTickLabel',{});
colormap(cbrewer('seq','Reds',20)); cb2=colorbar; freezeColors; cbfreeze(cb2);

h7=subplot(4,3,7);
contourf(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(cAllnorm0(50,:,:,31)),0.46:0.005:0.56)%,0.1:0.02:0.57)
caxis([0.46 0.56]); 
set(gca,'fontsize',12)
set(h7,'XTickLabel',{});
xlabel('(g)','fontsize',14)
ylabel({'Dye, t=299';'z'},'fontsize',14);
colormap(cbrewer('seq','Blues',20)); freezeColors;
h8=subplot(4,3,8);
contourf(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(cAllnorm(50,:,:,31)),0.46:0.005:0.56)%,0.1:0.02:0.57)
caxis([0.46 0.56]); 
set(gca,'fontsize',12)
colormap(cbrewer('seq','Blues',20)); freezeColors;
set(h8,'XTickLabel',{});
xlabel('(h)','fontsize',14)
set(h8,'YTickLabel',{});
h9=subplot(4,3,9);
contourf(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(cAllnorm2(50,:,:,31)),0.46:0.005:0.56)%,0.1:0.02:0.57)
caxis([0.46 0.56]); 
set(h9,'XTickLabel',{});
xlabel('(i)','fontsize',14)
set(h9,'YTickLabel',{});
set(gca,'fontsize',12)
colormap(cbrewer('seq','Blues',20)); cb3=colorbar; freezeColors; cbfreeze(cb3);

h10=subplot(4,3,10);
pcolor(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(keffmap0(50,:,:,31)))%,0.1:0.02:0.57)
shading 'flat'
caxis([0 0.016]); 
set(gca,'fontsize',12)
ylabel({'K_{eff}, t=299';'z'},'fontsize',14);
xlabel({'x';'(j)'},'fontsize',14); 
colormap(cbrewer('seq','Reds',20)); freezeColors;
h11=subplot(4,3,11);
pcolor(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(keffmap(50,:,:,31)))%,0.1:0.02:0.57)
shading 'flat'
caxis([0 0.016]); 
set(h11,'YTickLabel',{});
xlabel({'x';'(k)'},'fontsize',14); 
set(gca,'fontsize',12)
colormap(cbrewer('seq','Reds',20)); freezeColors;
h12=subplot(4,3,12);
pcolor(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(keffmap2(50,:,:,31)))%,0.1:0.02:0.57)
shading 'flat'
caxis([0 0.016]); 
xlabel({'x';'(l)'},'fontsize',14); 
set(h12,'YTickLabel',{});
set(gca,'fontsize',12)
colormap(cbrewer('seq','Reds',20)); cb4=colorbar; freezeColors; cbfreeze(cb4);
%% sept 2018 paper figure 8, full cylinder
load('E002compareLowDiffV2.mat','keffInt0','keffInt','time','chis0','chis','D_val')
figure; 
h1=subplot(2,3,1);
 plot(time,chis,'LineWidth',2); hold all; plot(time,chis0,'--','LineWidth',2)
ylabel('\chi^2','fontsize',14); xlabel('(a)','fontsize',14); title('E=0.02, k=10^{-6}','fontsize',16)
set(h1,'fontsize',14)
axis tight
set(h1,'XTickLabel',{});
h2=subplot(2,3,4);
plot(time,keffInt,'LineWidth',2); hold all; plot(time,keffInt0,'--','LineWidth',2); plot(time,ones(size(time))*(pi^6)*D_val./8,'k--')
xlabel({'time';'(d)'},'fontsize',14); ylabel('\int K_{eff} dV','fontsize',14);
set(h2,'fontsize',14); axis tight

load('E002compareV2.mat','keffInt0','keffInt','time','chis0','chis','D_val','keff0')
k0i=find(keff0>0,1,'first'); k0=keff0(k0i)*pi; keffInt0(1)=k0; keffInt(1)=k0;
h4=subplot(2,3,2);
plot(time,chis,'LineWidth',2); hold all; plot(time,chis0,'--','LineWidth',2)
set(h4,'XTickLabel',{}); 
xlabel('(b)','fontsize',14); title('E=0.02, k=10^{-4}','fontsize',16)
set(h4,'fontsize',14); axis tight
h3=subplot(2,3,5);
plot(time,keffInt,'LineWidth',2); hold all; plot(time,keffInt0,'--','LineWidth',2); plot(time,ones(size(time))*(pi^6)*D_val./8,'k--')
xlabel({'time','(e)'},'fontsize',14); 
set(h3,'fontsize',14); axis tight

load('E0125compare3v2.mat','keffInt0','keffInt','keffInt2','time','chis0','chis','chis2','D_val')
h5=subplot(2,3,6);
plot(time,keffInt,'LineWidth',2); hold all; plot(time,keffInt0,'--','LineWidth',2); plot(time,keffInt2,'LineWidth',2); 
plot(time,ones(size(time))*(pi^6)*D_val./8,'k--'); 
set(gca,'fontsize',14); axis tight
xlabel({'time';'(f)'},'fontsize',14)
legend('x0=-0.02','x0=0','x0=-0.16','analytic');

h6=subplot(2,3,3);
 plot(time,chis,'LineWidth',2); hold all; plot(time,chis0,'--','LineWidth',2); plot(time,chis2,'LineWidth',2)
xlabel('(c)','fontsize',14);
set(gca,'fontsize',14); axis tight
legend('x0=-0.02','x0=0','x0=-0.16'); 
set(h6,'XTickLabel',{}); 
title('E=0.125, k=10^{-4}','fontsize',16)
%% E=002 high-diff full cylinder and cross-sections
load('E002compareV2.mat')

figure; contourf(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(cAllnorm(50,:,:,5)),0.2:0.1:1)
caxis([0.2 1]); xlabel('x','fontsize',14); ylabel('z','fontsize',14); title('C, E=0.02, x_0=-0.02, k=10^{-4}, t=39','fontsize',14); colorbar
set(gca,'fontsize',12)
%contour 0.3:0.1:0.9
colormap(cbrewer('seq','Reds',8))

figure; contourf(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(cAllnorm(50,:,:,1)),0:0.1:1)
caxis([0 1]); xlabel('x','fontsize',14); ylabel('z','fontsize',14); title('C, E=0.02, x_0=-0.02, k=10^{-4}, t=39','fontsize',14); colorbar
set(gca,'fontsize',12)
%contour 0.3:0.1:0.9
colormap(cbrewer('seq','Reds',10))

figure; contourf(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(cAllnorm0(50,:,:,5)),0.3:0.1:1)
caxis([0.2 1]); xlabel('x','fontsize',14); ylabel('z','fontsize',14); title('C, E=0.02, x_0=0, k=10^{-4}, t=39','fontsize',14); colorbar
set(gca,'fontsize',12)
colormap(cbrewer('seq','Reds',8))

figure; contourf(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(cAllnorm(50,:,:,31)),0.455:0.005:0.545)
caxis([0.455 0.545]); 
xlabel('x','fontsize',14); ylabel('z','fontsize',14); title('C, E=0.02, x_0=-0.02, k=10^{-4}, t=299','fontsize',14); colorbar
set(gca,'fontsize',12)
colormap(cbrewer('seq','Reds',18))

figure; contourf(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(cAllnorm0(50,:,:,31)),0.455:0.005:0.545)
caxis([0.455 0.545]); 
xlabel('x','fontsize',14); ylabel('z','fontsize',14); title('C, E=0.02, x_0=0, k=10^{-4}, t=299','fontsize',14); colorbar
set(gca,'fontsize',12)
colormap(cbrewer('seq','Reds',18))

figure; pcolor(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(keffmap(50,:,:,5)))%,0:0.02:0.12)
caxis([0 0.12]); shading 'flat'
xlabel('x','fontsize',14); ylabel('z','fontsize',14); title('K, E=0.02, x_0=-0.02, k=10^{-4}, t=39','fontsize',14); colorbar
set(gca,'fontsize',12)

figure; pcolor(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(keffmap0(50,:,:,5)))%,0:0.02:0.12)
caxis([0 0.12]); shading 'flat'
xlabel('x','fontsize',14); ylabel('z','fontsize',14); title('K, E=0.02, x_0=0, k=10^{-4}, t=39','fontsize',14); colorbar
set(gca,'fontsize',12)

figure; pcolor(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(keffmap(50,:,:,31)))%,0:0.002:0.018)
caxis([0 0.018]); shading 'flat'
xlabel('x','fontsize',14); ylabel('z','fontsize',14); title('K, E=0.02, x_0=-0.02, k=10^{-4}, t=299','fontsize',14); colorbar
set(gca,'fontsize',12)

figure; pcolor(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(keffmap0(50,:,:,31)))%,0:0.002:0.018)
caxis([0 0.018]); shading 'flat'
xlabel('x','fontsize',14); ylabel('z','fontsize',14); title('K, E=0.02, x_0=0, k=10^{-4}, t=299','fontsize',14); colorbar
set(gca,'fontsize',12)

k0i=find(keff0>0,1,'first'); k0=keff0(k0i)*pi; keffInt0(1)=k0; keffInt(1)=k0;

figure; plot(time,keffInt,'LineWidth',2); hold all; plot(time,keffInt0,'--','LineWidth',2); plot(time,ones(size(time))*(pi^6)*D_val./8,'k--')
xlabel('time','fontsize',14); ylabel('\int K_{eff} dV','fontsize',14); legend('x0=-0.02','x0=0','t\rightarrow\infty'); title('E=0.02, k=10^{-4}','fontsize',14)
set(gca,'fontsize',12); axis tight

% figure; plot(time,keffInt/pi,'LineWidth',2); hold all; plot(time,keffInt0/pi,'LineWidth',2)
% xlabel('time','fontsize',14); ylabel('mean K_{eff}=\int K_{eff} dV/V','fontsize',14); legend('x0=-0.02','x0=0'); title('E=0.02, k=10^{-4}','fontsize',14)
% set(gca,'fontsize',12); axis tight

figure; plot(time,chis,'LineWidth',2); hold all; plot(time,chis0,'--','LineWidth',2)
xlabel('time','fontsize',14); ylabel('\chi^2','fontsize',14); legend('x0=-0.02','x0=0'); title('E=0.02, k=10^{-4}','fontsize',14)
set(gca,'fontsize',12); axis tight

% keffmap1=keffmap;
% keffmap01=keffmap0;
% for i=1:nt
%    holdvar=keffmap1(:,:,:,i); holdvar(inside==0)=nan; keffmap1(:,:,:,i)=holdvar;
%    holdvar=keffmap01(:,:,:,i); holdvar(inside==0)=nan; keffmap01(:,:,:,i)=holdvar;
% end
% ratio1=keffmap1./keffmap01; clear keffmap1 keffmap01
% ratio1(isinf(ratio1))=NaN;
% rKmean=squeeze(nanmean(nanmean(nanmean(ratio1))));
% figure; plot(time,rKmean,'LineWidth',2); title('Mean ratio K_{chaos}/K_{reg}, E=0.02, k=10^{-4}','fontsize',14)
% xlabel('time','fontsize',14); ylabel('R_K','fontsize',14); set(gca,'fontsize',12)
%% E002 high diff regions, updated Sept 2018 to make rcFigure11
load('E002compareV2.mat')

x1=squeeze(x(50,:,:)); z1=squeeze(z(50,:,:));
y2=squeeze(y(:,50,:)); z2=squeeze(z(:,50,:));

polyx1=[-0.55 -0.55 -0.45 -0.45 -0.4 -0.4 -0.35 -0.35 -0.4 -0.4 -0.45 -0.45 -0.55];
polyz1=[0.35 0.45 0.45 0.6 0.6 0.65 0.65 0.45 0.45 0.4 0.4 0.35 0.35];

inisland=double(inpolygon(x1,z1,polyx1,polyz1));
polyx2=[0.55 0.7 0.7 0.75 0.75 0.8 0.8 0.75 0.75 0.7 0.7 0.65 0.65 0.55 0.55];
polyz2=[0.85 0.85 0.8 0.8 0.75 0.75 0.65 0.65 0.7 0.7 0.75 0.75 0.8 0.8 0.85];
inisland=inisland+inpolygon(x1,z1,polyx2,polyz2);

inchaos2=double(inpolygon(x1,z1,-polyx1,polyz1));
inchaos2=inchaos2+double(inpolygon(x1,z1,-polyx2,polyz2));

%clear kisland2* kchaos2* cisland2* cchaos2*
for i=1:41
holdvar1=squeeze(cAllnorm(50,:,:,i));
holdvar2=holdvar1(inisland==1); cisland2(i,:)=holdvar2;
holdvar1=squeeze(cAllnorm0(50,:,:,i));
holdvar2=holdvar1(inisland==1); cisland20(i,:)=holdvar2;

holdvar1=squeeze(keffmap(50,:,:,i));
holdvar2=holdvar1(inisland==1); kisland2(i,:)=holdvar2;
holdvar1=squeeze(keffmap0(50,:,:,i));
holdvar2=holdvar1(inisland==1); kisland20(i,:)=holdvar2;

holdvar1=squeeze(cAllnorm(50,:,:,i));
holdvar2=holdvar1(inchaos2==1); cchaos2(i,:)=holdvar2;
holdvar1=squeeze(cAllnorm0(50,:,:,i));
holdvar2=holdvar1(inchaos2==1); cchaos20(i,:)=holdvar2;

holdvar1=squeeze(keffmap(50,:,:,i));
holdvar2=holdvar1(inchaos2==1); kchaos2(i,:)=holdvar2;
holdvar1=squeeze(keffmap0(50,:,:,i));
holdvar2=holdvar1(inchaos2==1); kchaos20(i,:)=holdvar2;
end

polyy1=[-.65 -0.65 -0.45 -0.45 -0.35 -0.35 -0.45 -0.45 -0.5 -0.5 -0.65];
polyzy1=[0.75 0.85 0.85 0.8 0.8 0.5 0.5 0.7 0.7 0.75 0.75];
polyy2=[0.45 0.45 0.65 0.65 0.7 0.7 0.8 0.8 0.72 0.72 0.68 0.68 0.45];
polyzy2=[0.35 0.45 0.45 0.5 0.5 0.7 0.7 0.5 0.5 0.4 0.4 0.35 0.35];

inislandy=double(inpolygon(y2,z2,polyy1,polyzy1));
inislandy=inislandy+double(inpolygon(y2,z2,polyy2,polyzy2));
inchaos2y=double(inpolygon(y2,z2,-polyy1,polyzy1));
inchaos2y=inchaos2y+double(inpolygon(y2,z2,-polyy2,polyzy2));

for i=1:41
holdvar1=squeeze(cAllnorm(:,50,:,i));
holdvar2=holdvar1(inislandy==1); cislandy(i,:)=holdvar2;
holdvar1=squeeze(cAllnorm0(:,50,:,i));
holdvar2=holdvar1(inislandy==1); cislandy0(i,:)=holdvar2;

holdvar1=squeeze(keffmap(:,50,:,i));
holdvar2=holdvar1(inislandy==1); kislandy(i,:)=holdvar2;
holdvar1=squeeze(keffmap0(:,50,:,i));
holdvar2=holdvar1(inislandy==1); kislandy0(i,:)=holdvar2;

holdvar1=squeeze(cAllnorm(:,50,:,i));
holdvar2=holdvar1(inchaos2y==1); cchaos2y(i,:)=holdvar2;
holdvar1=squeeze(cAllnorm0(:,50,:,i));
holdvar2=holdvar1(inchaos2y==1); cchaos2y0(i,:)=holdvar2;

holdvar1=squeeze(keffmap(:,50,:,i));
holdvar2=holdvar1(inchaos2y==1); kchaos2y(i,:)=holdvar2;
holdvar1=squeeze(keffmap0(:,50,:,i));
holdvar2=holdvar1(inchaos2y==1); kchaos2y0(i,:)=holdvar2;
end

kisland=cat(2,kisland2,kislandy);
kisland0=cat(2,kisland20,kislandy0);
kchaos=cat(2,kchaos2,kchaos2y);
kchaos0=cat(2,kchaos20,kchaos2y0);

risland=kisland./kisland0; risland(isinf(risland))=nan;
rchaos=kchaos./kchaos0; rchaos(isinf(rchaos))=nan;

figure; h1=subplot(2,2,1);
 plot(-polyx2,polyz2,'r','LineWidth',2)
hold on; plot(-polyx1,polyz1,'r','LineWidth',2)
plot(polyx2,polyz2,'b','LineWidth',2)
plot(polyx1,polyz1,'b','LineWidth',2)
scatter(poincareX,poincareZ,1,'k')
xlabel('x','fontsize',14); ylabel({'Poincare';'z'},'fontsize',14); 
set(h1,'fontsize',14)
h2=subplot(2,2,2);
plot(-polyy2,polyzy2,'r','LineWidth',2)
hold on; plot(-polyy1,polyzy1,'r','LineWidth',2)
plot(polyy2,polyzy2,'b','LineWidth',2)
plot(polyy1,polyzy1,'b','LineWidth',2)
scatter(poincareY,poincareZy,1,'k')
xlabel('y','fontsize',14); set(h2,'fontsize',12)
set(h2,'YTickLabel',{});
h3=subplot(2,2,3);
plot(time,mean(kchaos,2),'r','LineWidth',2)
hold on
plot(time,mean(kisland,2),'b','LineWidth',2)
plot(time,mean(kchaos0,2),'Color',[0 0.7 0],'LineWidth',2)
legend('resonant','island','no-chaos match')
ylabel('K_{eff}','fontsize',14); xlabel({'time';'k=10^{-4}'},'fontsize',14); 
set(h3,'fontsize',14); axis tight

load('E002compareLowDiffV2.mat')
load('E002compareWorkingSep23.mat', 'in*')

kxz=squeeze(keffmap(50,:,:,:));
kxz0=squeeze(keffmap0(50,:,:,:));
kyz=squeeze(keffmap(:,50,:,:));
kyz0=squeeze(keffmap0(:,50,:,:));

inislandy=squeeze(inislandy);
inchaos2y=squeeze(inchaos2y);

for i=1:41
    holdvar1=kxz(:,:,i);
    holdvar2=holdvar1(inisland==1); kislandx(i,:)=holdvar2;
    holdvar2=holdvar1(inchaos2==1); kchaosx(i,:)=holdvar2;
    %holdvar2=holdvar1(incenterxz==1); kcenterx(i,:)=holdvar2;
    %holdvar2=holdvar1(inouterxz==1); kouterx(i,:)=holdvar2;
    
    holdvar1=squeeze(kxz0(:,:,i));
    holdvar2=holdvar1(inisland==1); kislandx0(i,:)=holdvar2;
    holdvar2=holdvar1(inchaos2==1); kchaosx0(i,:)=holdvar2;
    %holdvar2=holdvar1(incenterxz==1); kcenterx0(i,:)=holdvar2;
    %holdvar2=holdvar1(inouterxz==1); kouterx0(i,:)=holdvar2;
    
    holdvar1=kyz(:,:,i);
    holdvar2=holdvar1(inislandy==1); kislandy(i,:)=holdvar2;
    holdvar2=holdvar1(inchaos2y==1); kchaosy(i,:)=holdvar2;
    %holdvar2=holdvar1(incenteryz==1); kcentery(i,:)=holdvar2;
    %holdvar2=holdvar1(inouteryz==1); koutery(i,:)=holdvar2;
    
    holdvar1=squeeze(kyz0(:,:,i));
    holdvar2=holdvar1(inislandy==1); kislandy0(i,:)=holdvar2;
    holdvar2=holdvar1(inchaos2y==1); kchaosy0(i,:)=holdvar2;
    %holdvar2=holdvar1(incenteryz==1); kcentery0(i,:)=holdvar2;
    %holdvar2=holdvar1(inouteryz==1); koutery0(i,:)=holdvar2;
end

%kouter=cat(2,kouterx,koutery);
%kouter0=cat(2,kouterx0,koutery0);
%kcenter=cat(2,kcenterx,kcentery);
%kcenter0=cat(2,kcenterx0,kcentery0);
kisland=cat(2,kislandx,kislandy);
kisland0=cat(2,kislandx0,kislandy0);
kchaos=cat(2,kchaosx,kchaosy);
kchaos0=cat(2,kchaosx0,kchaosy0);

h4=subplot(2,2,4);
plot(time,mean(kchaos,2),'r','LineWidth',2)
hold on
plot(time,mean(kisland,2),'b','LineWidth',2)
plot(time,mean(kchaos0,2),'Color',[0 0.7 0],'LineWidth',2)
xlabel({'time';'k=10^{-6}'},'fontsize',14); 
set(h4,'fontsize',14); axis tight



%% E=0.02 low-diff full cylinder and cross-sections
%load('E002compareLowDiffV2.mat')
figure; contourf(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(cAllnorm(50,:,:,5)),0.2:0.1:1)
caxis([0.2 1]); xlabel('x','fontsize',14); ylabel('z','fontsize',14); title('C, E=0.02, x_0=-0.02, k=10^{-6}, t=39','fontsize',14); colorbar
set(gca,'fontsize',12)
colormap(cbrewer('seq','Reds',8))

figure; contourf(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(cAllnorm0(50,:,:,5)),0:0.1:1)
caxis([0 1]); xlabel('x','fontsize',14); ylabel('z','fontsize',14); title('C, E=0.02, x_0=0, k=10^{-6}, t=39','fontsize',14); colorbar
set(gca,'fontsize',12)
colormap(cbrewer('seq','Reds',10))
%
figure; contourf(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(cAllnorm(50,:,:,31)),0.3:0.05:0.7)
%caxis([0 1]); 
xlabel('x','fontsize',14); ylabel('z','fontsize',14); title('C, E=0.02, x_0=-0.02, k=10^{-6}, t=299','fontsize',14); colorbar
set(gca,'fontsize',11)
colormap(cbrewer('seq','Reds',8))

figure; contourf(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(cAllnorm0(50,:,:,31)),0.3:0.05:0.7)
%caxis([0.1 0.57]); 
xlabel('x','fontsize',14); ylabel('z','fontsize',14); title('C, E=0.02, x_0=0, k=10^{-6}, t=299','fontsize',14); colorbar
set(gca,'fontsize',11)
colormap(cbrewer('seq','Reds',8))

%
figure; pcolor(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(keffmap(50,:,:,5)))%,0:0.05:0.27)
caxis([0 0.25]); shading 'flat'
xlabel('x','fontsize',14); ylabel('z','fontsize',14); title('K, E=0.02, x_0=-0.02, k=10^{-6}, t=39','fontsize',14); colorbar
set(gca,'fontsize',12)

figure; pcolor(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(keffmap0(50,:,:,5)))%,0:0.05:0.27)
caxis([0 0.25]); shading 'flat'
xlabel('x','fontsize',14); ylabel('z','fontsize',14); title('K, E=0.02, x_0=0, k=10^{-6}, t=39','fontsize',14); colorbar
set(gca,'fontsize',12)

figure; pcolor(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(keffmap(50,:,:,31)))%,[0 0.005 0.01:0.01:0.1])
caxis([0 0.1]); shading 'flat'
xlabel('x','fontsize',14); ylabel('z','fontsize',14); title('K, E=0.02, x_0=-0.02, k=10^{-6}, t=299','fontsize',14); colorbar
set(gca,'fontsize',12)

figure; pcolor(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(keffmap0(50,:,:,31)))%,[0 0.005 0.01:0.01:0.1])
caxis([0 0.1]); shading 'flat'
xlabel('x','fontsize',14); ylabel('z','fontsize',14); title('K, E=0.02, x_0=0, k=10^{-6}, t=299','fontsize',14); colorbar
set(gca,'fontsize',12)

figure; plot(time,keffInt,'LineWidth',2); hold all; plot(time,keffInt0,'--','LineWidth',2); plot(time,ones(size(time))*(pi^6)*D_val./8,'k--')
xlabel('time','fontsize',14); ylabel('\int K_{eff} dV','fontsize',14); legend('x0=-0.02','x0=0','t\rightarrow\infty'); title('E=0.02, k=10^{-6}','fontsize',14)
set(gca,'fontsize',12); axis tight

figure; plot(time,chis,'LineWidth',2); hold all; plot(time,chis0,'--','LineWidth',2)
xlabel('time','fontsize',14); ylabel('\chi^2','fontsize',14); legend('x0=-0.02','x0=0'); title('E=0.02, k=10^{-6}','fontsize',14)
set(gca,'fontsize',12); axis tight

keffmap1=keffmap;
keffmap01=keffmap0;
for i=1:nt
   holdvar=keffmap1(:,:,:,i); holdvar(inside==0)=nan; keffmap1(:,:,:,i)=holdvar;
   holdvar=keffmap01(:,:,:,i); holdvar(inside==0)=nan; keffmap01(:,:,:,i)=holdvar;
end
ratio1=keffmap1./keffmap01; clear keffmap1 keffmap01
ratio1(isinf(ratio1))=NaN;
rKmean=squeeze(nanmean(nanmean(nanmean(ratio1))));
figure; semilogy(time,rKmean,'LineWidth',2); title('Mean ratio K_{chaos}/K_{reg}, E=0.02, k=10^{-6}','fontsize',14)
xlabel('time','fontsize',14); ylabel('R_K','fontsize',14); set(gca,'fontsize',12); axis tight
%% E=0.02 low-diff regions 
 load('E002compareWorkingSep23.mat', 'in*')

kxz=squeeze(keffmap(50,:,:,:));
kxz0=squeeze(keffmap0(50,:,:,:));
kyz=squeeze(keffmap(:,50,:,:));
kyz0=squeeze(keffmap0(:,50,:,:));

inislandy=squeeze(inislandy);
inchaos2y=squeeze(inchaos2y);

for i=1:41
    holdvar1=kxz(:,:,i);
    holdvar2=holdvar1(inisland==1); kislandx(i,:)=holdvar2;
    holdvar2=holdvar1(inchaos2==1); kchaosx(i,:)=holdvar2;
    holdvar2=holdvar1(incenterxz==1); kcenterx(i,:)=holdvar2;
    holdvar2=holdvar1(inouterxz==1); kouterx(i,:)=holdvar2;
    
    holdvar1=squeeze(kxz0(:,:,i));
    holdvar2=holdvar1(inisland==1); kislandx0(i,:)=holdvar2;
    holdvar2=holdvar1(inchaos2==1); kchaosx0(i,:)=holdvar2;
    holdvar2=holdvar1(incenterxz==1); kcenterx0(i,:)=holdvar2;
    holdvar2=holdvar1(inouterxz==1); kouterx0(i,:)=holdvar2;
    
    holdvar1=kyz(:,:,i);
    holdvar2=holdvar1(inislandy==1); kislandy(i,:)=holdvar2;
    holdvar2=holdvar1(inchaos2y==1); kchaosy(i,:)=holdvar2;
    holdvar2=holdvar1(incenteryz==1); kcentery(i,:)=holdvar2;
    holdvar2=holdvar1(inouteryz==1); koutery(i,:)=holdvar2;
    
    holdvar1=squeeze(kyz0(:,:,i));
    holdvar2=holdvar1(inislandy==1); kislandy0(i,:)=holdvar2;
    holdvar2=holdvar1(inchaos2y==1); kchaosy0(i,:)=holdvar2;
    holdvar2=holdvar1(incenteryz==1); kcentery0(i,:)=holdvar2;
    holdvar2=holdvar1(inouteryz==1); koutery0(i,:)=holdvar2;
end

kouter=cat(2,kouterx,koutery);
kouter0=cat(2,kouterx0,koutery0);
kcenter=cat(2,kcenterx,kcentery);
kcenter0=cat(2,kcenterx0,kcentery0);
kisland=cat(2,kislandx,kislandy);
kisland0=cat(2,kislandx0,kislandy0);
kchaos=cat(2,kchaosx,kchaosy);
kchaos0=cat(2,kchaosx0,kchaosy0);

router=kouter./kouter0; router(isinf(router))=nan;

figure; plot(time,nanmean(kchaos./kchaos0,2),time,nanmean(kisland./kisland0,2),'LineWidth',2)%,time,2*ones(41,1),'k')
xlabel('time','fontsize',14); ylabel('R_K','fontsize',14); title('Mean R_{K} in different regions','fontsize',14)
legend('resonant','island'); set(gca,'fontsize',12); axis tight

figure; plot(time,nanmean(kchaos,2),time,nanmean(kisland,2),'LineWidth',2)
hold on; plot(time,nanmean(kchaos0,2),':',time,nanmean(kisland0,2),':','LineWidth',2)
xlabel('time','fontsize',14); ylabel('K_{eff}','fontsize',14); title('Mean K_{eff} in different regions, E=0.02, k=10^{-6}','fontsize',14)
legend('resonant','island','no-chaos match'); set(gca,'fontsize',12); axis tight


%% E=0.125 full cylinder and cross-sections
load('E0125compare3v2.mat')
figure; contourf(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(cAllnorm(50,:,:,5)),0.2:0.1:1)
caxis([0.2 1]); xlabel('x','fontsize',14); ylabel('z','fontsize',14); title('C, E=0.125, x_0=-0.02, k=10^{-4}, t=39','fontsize',14); colorbar
set(gca,'fontsize',12)
colormap(cbrewer('seq','Reds',8))

figure; contourf(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(cAllnorm0(50,:,:,5)))
caxis([0.2 1]); xlabel('x','fontsize',14); ylabel('z','fontsize',14); title('C, E=0.125, x_0=0, k=10^{-4}, t=39','fontsize',14); colorbar
set(gca,'fontsize',12)
colormap(cbrewer('seq','Reds',8))

figure; contourf(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(cAllnorm2(50,:,:,5)))
caxis([0.2 1]); xlabel('x','fontsize',14); ylabel('z','fontsize',14); title('C, E=0.125, x_0=-0.16, k=10^{-4}, t=39','fontsize',14); colorbar
set(gca,'fontsize',12)
colormap(cbrewer('seq','Reds',8))

figure; contourf(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(cAllnorm(50,:,:,31)),0.46:0.005:0.56)%,0.1:0.02:0.57)
caxis([0.46 0.56]); 
xlabel('x','fontsize',14); ylabel('z','fontsize',14); title('C, E=0.125, x_0=-0.02, k=10^{-4}, t=299','fontsize',14); colorbar
set(gca,'fontsize',12)
colormap(cbrewer('seq','Reds',20))

figure; contourf(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(cAllnorm0(50,:,:,31)),0.46:0.005:0.56)%,0.15:0.02:0.57)
caxis([0.46 0.56]); 
xlabel('x','fontsize',14); ylabel('z','fontsize',14); title('C, E=0.125, x_0=0, k=10^{-4}, t=299','fontsize',14); colorbar
set(gca,'fontsize',12)
colormap(cbrewer('seq','Reds',20))

figure; contourf(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(cAllnorm2(50,:,:,31)),0.46:0.005:0.56)%,0.15:0.02:0.57)
caxis([0.46 0.56]); 
xlabel('x','fontsize',14); ylabel('z','fontsize',14); title('C, E=0.125, x_0=-0.16, k=10^{-4}, t=299','fontsize',14); colorbar
set(gca,'fontsize',12)
colormap(cbrewer('seq','Reds',20))
%%
figure; pcolor(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(keffmap(50,:,:,5)))%,0:0.002:0.016)%,0:0.05:0.27)
caxis([0 0.015]); shading 'flat'
xlabel('x','fontsize',14); ylabel('z','fontsize',14); title('K, E=0.125, x_0=-0.02, k=10^{-4}, t=39','fontsize',14); colorbar
set(gca,'fontsize',12)

figure; pcolor(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(keffmap0(50,:,:,5)))%,0:0.003:0.016)%,0:0.05:0.27)
caxis([0 0.015]); shading 'flat'
xlabel('x','fontsize',14); ylabel('z','fontsize',14); title('K, E=0.125, x_0=0, k=10^{-4}, t=39','fontsize',14); colorbar
set(gca,'fontsize',12)

figure; pcolor(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(keffmap2(50,:,:,5)))%,0:0.003:0.016)
caxis([0 0.015]); shading 'flat'
xlabel('x','fontsize',14); ylabel('z','fontsize',14); title('K, E=0.125, x_0=-0.16, k=10^{-4}, t=39','fontsize',14); colorbar
set(gca,'fontsize',12)

figure; pcolor(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(keffmap(50,:,:,31)))%,[0 0.0001 0.00025 0.001:0.002:0.014])%,[0 0.005 0.01:0.01:0.1])
caxis([0 0.016]); shading 'flat'
xlabel('x','fontsize',14); ylabel('z','fontsize',14); title('K, E=0.125, x_0=-0.02, k=10^{-4}, t=299','fontsize',14); colorbar
set(gca,'fontsize',12)

figure; pcolor(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(keffmap0(50,:,:,31)))%,[0 0.0001 0.00025 0.001:0.002:0.014])%,[0 0.005 0.01:0.01:0.1])
caxis([0 0.016]); shading 'flat'
xlabel('x','fontsize',14); ylabel('z','fontsize',14); title('K, E=0.125, x_0=0, k=10^{-4}, t=299','fontsize',14); colorbar
set(gca,'fontsize',12)

figure; pcolor(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(keffmap2(50,:,:,31)))%,[0 0.0001 0.00025 0.001:0.002:0.014])%%caxis([0 0.025]); 
caxis([0 0.016]); shading 'flat'
xlabel('x','fontsize',14); ylabel('z','fontsize',14); title('K, E=0.125, x_0=-0.16, k=10^{-4}, t=299','fontsize',14); colorbar
set(gca,'fontsize',12)

figure; plot(time,keffInt,'LineWidth',2); hold all; plot(time,keffInt0,'--','LineWidth',2); plot(time,keffInt2,'LineWidth',2); 
plot(time,ones(size(time))*(pi^6)*D_val./8,'k--'); plot(time,ones(size(time))*8*(pi^3)*D_val,'r--')
legend('x0=-0.02','x0=0','x0=-0.16','t\rightarrow\infty, ct','t\rightarrow\infty, st');
xlabel('time','fontsize',14); ylabel('\int K_{eff} dV','fontsize',14); title('E=0.125, k=10^{-4}','fontsize',14)
set(gca,'fontsize',12); axis tight

figure; plot(time,keffInt/pi,'LineWidth',2); hold all; plot(time,keffInt0/pi,'LineWidth',2); plot(time,keffInt2/pi,'LineWidth',2)
xlabel('time','fontsize',14); ylabel('mean K_{eff}','fontsize',14); legend('x0=-0.02','x0=0','x0=-0.16'); title('E=0.125, k=10^{-4}','fontsize',14)
set(gca,'fontsize',12); axis tight

figure; plot(time,chis,'LineWidth',2); hold all; plot(time,chis0,'--','LineWidth',2); plot(time,chis2,'LineWidth',2)
xlabel('time','fontsize',14); ylabel('\chi^2','fontsize',14); legend('x0=-0.02','x0=0','x0=-0.16'); title('E=0.125, k=10^{-4}','fontsize',14)
set(gca,'fontsize',12); axis tight

keffmap1=keffmap;
keffmap01=keffmap0;
keffmap21=keffmap2;
for i=1:nt
   holdvar=keffmap1(:,:,:,i); holdvar(inside==0)=nan; keffmap1(:,:,:,i)=holdvar;
   holdvar=keffmap01(:,:,:,i); holdvar(inside==0)=nan; keffmap01(:,:,:,i)=holdvar;
   holdvar=keffmap21(:,:,:,i); holdvar(inside==0)=nan; keffmap21(:,:,:,i)=holdvar;
end
ratio1=keffmap1./keffmap01; ratio2=keffmap21./keffmap01; clear keffmap1 keffmap01 keffmap21
ratio1(isinf(ratio1))=NaN; ratio2(isinf(ratio2))=NaN;
rKmean=squeeze(nanmean(nanmean(nanmean(ratio1)))); rKmean2=squeeze(nanmean(nanmean(nanmean(ratio2))));
figure; semilogy(time,rKmean,'LineWidth',2); hold all; semilogy(time,rKmean2,'LineWidth',2);
title('Mean ratio K_{chaos}/K_{reg}, E=0.02, k=10^{-4}','fontsize',14)
xlabel('time','fontsize',14); ylabel('R_K','fontsize',14); set(gca,'fontsize',12); axis tight
legend('x0=-0.02','x0=-0.16')

%% E=0.125 regions x0=-0.16

x11=[-0.52 -0.43 -0.34 -0.4 -0.53];
zx11=[0.52 0.51 0.39 0.35 0.43];
x12=[0.46 0.54 0.61 0.55 0.45];
zx12=[0.81 0.825 0.805 0.75 0.76];
zx21=[0.68 0.67 0.3 0.22 0.25 0.35 0.45 0.65];
x22=[0.45 0.65 0.74 0.74 0.56 0.4 0.32 0.32];
zx22=[0.88 0.88 0.84 0.8 0.68 0.68 0.76 0.8];
x21=[-0.55 -0.48 -0.1 -0.2 -0.4 -0.6 -0.65 -0.65];

y11=[-0.55 -0.46 -0.43 -0.52 -0.6];
zy11=[0.72 0.7 0.62 0.59 0.67];
y12=[0.57 0.65 0.63 0.54 0.53];
zy12=[0.68 0.66 0.5 0.48 0.62];
y21=[-0.65 -0.45 -0.35 -0.3 -0.35 -0.5 -0.7 -0.73];
zy21=[0.82 0.82 0.75 0.5 0.4 0.4 0.6 0.7];
y22=[0.45 0.7 0.75 0.75 0.55 0.35 0.3 0.4 0.4];
zy22=[0.8 0.8 0.7 0.5 0.2 0.2 0.3 0.45 0.7];

incenterx=double(inpolygon(squeeze(x(50,:,:)),squeeze(z(50,:,:)),x11,zx11))+double(inpolygon(squeeze(x(50,:,:)),squeeze(z(50,:,:)),x12,zx12));
incentery=double(inpolygon(squeeze(y(:,50,:)),squeeze(z(:,50,:)),y11,zy11))+double(inpolygon(squeeze(y(:,50,:)),squeeze(z(:,50,:)),y12,zy12));
inresy=double(inpolygon(squeeze(y(:,50,:)),squeeze(z(:,50,:)),y21,zy21))+double(inpolygon(squeeze(y(:,50,:)),squeeze(z(:,50,:)),y22,zy22))-incentery;
inresx=double(inpolygon(squeeze(x(50,:,:)),squeeze(z(50,:,:)),x21,zx21))+double(inpolygon(squeeze(x(50,:,:)),squeeze(z(50,:,:)),x22,zx22))-incenterx;
inouterx=ones(size(incenterx))-incenterx-inresx;
inoutery=ones(size(incentery))-incentery-inresy;

for i=1:41
    holdvar1=squeeze(keffmap2(50,:,:,i));
    holdvar2=holdvar1(incenterx==1); kcenx(i,:)=holdvar2;
    holdvar2=holdvar1(inresx==1); kresx(i,:)=holdvar2;
    holdvar2=holdvar1(inouterx==1); koutx(i,:)=holdvar2;
    
    holdvar1=squeeze(keffmap0(50,:,:,i));
    holdvar2=holdvar1(incenterx==1); kcenx0(i,:)=holdvar2;
    holdvar2=holdvar1(inresx==1); kresx0(i,:)=holdvar2;
    holdvar2=holdvar1(inouterx==1); koutx0(i,:)=holdvar2;

    holdvar1=squeeze(keffmap2(:,50,:,i));
    holdvar2=holdvar1(incentery==1); kceny(i,:)=holdvar2;
    holdvar2=holdvar1(inresy==1); kresy(i,:)=holdvar2;
    holdvar2=holdvar1(inoutery==1); kouty(i,:)=holdvar2;
    
    holdvar1=squeeze(keffmap0(:,50,:,i));
    holdvar2=holdvar1(incentery==1); kceny0(i,:)=holdvar2;
    holdvar2=holdvar1(inresy==1); kresy0(i,:)=holdvar2;
    holdvar2=holdvar1(inoutery==1); kouty0(i,:)=holdvar2;
    
end

kout=cat(2,koutx,kouty);
kout0=cat(2,koutx0,kouty0);
kcen=cat(2,kcenx,kceny);
kcen0=cat(2,kcenx0,kceny0);
kres=cat(2,kresx,kresy);
kres0=cat(2,kresx0,kresy0);

router=kout./kout0; router(isinf(router))=nan;
rres=kres./kres0; rres(isinf(rres))=nan;
rcenter=kcen./kcen0; rcenter(isinf(rcenter))=nan;

% figure; semilogy(time,nanmean(rres,2),time,nanmean(rcenter,2),time,nanmean(router,2))%,time,2*ones(41,1),'k')
% xlabel('time','fontsize',14); ylabel('R_K','fontsize',14); title('Mean R_{K} in different regions, E=0.125, x_0=-0.16','fontsize',14)
% legend('resonant','center','outer')
% 
% figure; plot(time,nanmean(kres,2),time,nanmean(kcen,2),time,nanmean(kout,2),'LineWidth',2)
% hold on; plot(time,nanmean(kres0,2),':',time,nanmean(kcen0,2),':',time,nanmean(kout0,2),':','LineWidth',2)
% xlabel('time','fontsize',14); ylabel('K_{eff}','fontsize',14); title('Mean K_{eff} in different regions, E=0.125, x_0=-0.16','fontsize',14)
% legend('resonant','center','outer','resonant regular','center regular','outer regular')
% 
% figure; plot(time,nanmean(kres,2)./nanmean(kcen,2),time,nanmean(kres0,2)./nanmean(kcen0,2),'LineWidth',2)
% legend('perturbed','unperturbed')
% axis([0 400 0 50])
% title('Mean K_{resonant}/K_{center}','fontsize',14)
% hold on; plot(time,2*ones(size(time)),'k'); axis tight

figure; plot(time,nanmean(kres,2),time,nanmean(kcen,2),'LineWidth',2)
hold on; plot(time,nanmean(kcen0,2),':','LineWidth',2)
xlabel('time','fontsize',14); ylabel('K_{eff}','fontsize',14); title('Mean K_{eff} in different regions, E=0.125, x_0=-0.16','fontsize',14)
legend('resonant','center','center regular'); set(gca,'fontsize',12); axis tight

%% E=0.125 regions x0=-0.02

load('E0125compareWorkingSep22.mat','in*')

kxz=squeeze(keffmap(50,:,:,:));
kxz0=squeeze(keffmap0(50,:,:,:));
kyz=squeeze(keffmap(:,50,:,:));
kyz0=squeeze(keffmap0(:,50,:,:));

for i=1:41
    holdvar1=kxz(:,:,i);
    holdvar2=holdvar1(inislxz==1); kislandx(i,:)=holdvar2;
    holdvar2=holdvar1(inresxz==1); kchaosx(i,:)=holdvar2;
    holdvar2=holdvar1(incenterxz==1); kcenterx(i,:)=holdvar2;
    holdvar2=holdvar1(inouterxz==1); kouterx(i,:)=holdvar2;
    
    holdvar1=squeeze(kxz0(:,:,i));
    holdvar2=holdvar1(inislxz==1); kislandx0(i,:)=holdvar2;
    holdvar2=holdvar1(inresxz==1); kchaosx0(i,:)=holdvar2;
    holdvar2=holdvar1(incenterxz==1); kcenterx0(i,:)=holdvar2;
    holdvar2=holdvar1(inouterxz==1); kouterx0(i,:)=holdvar2;
    
    holdvar1=kyz(:,:,i);
    holdvar2=holdvar1(inislyz==1); kislandy(i,:)=holdvar2;
    holdvar2=holdvar1(inresyz==1); kchaosy(i,:)=holdvar2;
    holdvar2=holdvar1(incenteryz==1); kcentery(i,:)=holdvar2;
    holdvar2=holdvar1(inouteryz==1); koutery(i,:)=holdvar2;
    
    holdvar1=squeeze(kyz0(:,:,i));
    holdvar2=holdvar1(inislyz==1); kislandy0(i,:)=holdvar2;
    holdvar2=holdvar1(inresyz==1); kchaosy0(i,:)=holdvar2;
    holdvar2=holdvar1(incenteryz==1); kcentery0(i,:)=holdvar2;
    holdvar2=holdvar1(inouteryz==1); koutery0(i,:)=holdvar2;
end

kouter=cat(2,kouterx,koutery);
kouter0=cat(2,kouterx0,koutery0);
kcenter=cat(2,kcenterx,kcentery);
kcenter0=cat(2,kcenterx0,kcentery0);
kisland=cat(2,kislandx,kislandy);
kisland0=cat(2,kislandx0,kislandy0);
kchaos=cat(2,kchaosx,kchaosy);
kchaos0=cat(2,kchaosx0,kchaosy0);

router=kouter./kouter0; router(isinf(router))=nan;

figure; plot(time,nanmean(kchaos./kchaos0,2),time,nanmean(kisland./kisland0,2),'LineWidth',2)%,time,2*ones(41,1),'k')
xlabel('time','fontsize',14); ylabel('R_K','fontsize',14); title('Mean R_{K} in different regions','fontsize',14)
legend('resonant','island'); set(gca,'fontsize',12); axis tight

figure; plot(time,nanmean(kchaos,2),time,nanmean(kisland,2),'LineWidth',2)
hold on; plot(time,nanmean(kchaos0,2),':',time,nanmean(kisland0,2),':','LineWidth',2)
xlabel('time','fontsize',14); ylabel('K_{eff}','fontsize',14); title('Mean K_{eff} in different regions','fontsize',14)
legend('resonant','island','resonant regular','island regular'); set(gca,'fontsize',12); axis tight

%% kefflines
load('E002compareV2.mat')
keffE02=keff; keff0E02=keff0;
load('E002compareLowDiffV2.mat')
keffE02low=keff; keff0E02low=keff0;
load('E0125compare3v2.mat')

%% early time, t=39, i=5
figure; plot(cvaldvdc,keffE02low(5,:),'linewidth',2); hold all;
plot(cvaldvdc,keff0E02low(5,:),'--')
plot(cvaldvdc,keffE02(5,:),'linewidth',2)
plot(cvaldvdc,keff0E02(5,:),'--')
%plot(cvaldvdc,keff(5,:),'linewidth',2)
%plot(cvaldvdc,keff2(5,:),'linewidth',2)
%plot(cvaldvdc,keff0(5,:),'--')
axis([0 1 0 0.5])
%legend('E=0.02,k=10^{-6},x_0=-0.02','E=0.02,k=10^{-6},x_0=0','E=0.02,k=10^{-4},x_0=-0.02','E=0.02,k=10^{-4},x_0=0','E=0.125,k=10^{-4},x_0=-0.02','E=0.125,k=10^{-4},x_0=-0.16','E=0.125,k=10^{-4},x_0=0')
set(gca,'fontsize',22)
%set(gca,'yaxis','log')
xlabel('dye concentration','fontsize',22)
ylabel('\kappa_{eff}','fontsize',22)
%% late time
max1=max(keffE02(31,:));
maxi1=find(keffE02==max1);
figure; plot(cvaldvdc,keffE02low(31,:),'linewidth',2); hold all;
plot(cvaldvdc,keff0E02low(31,:),'--')
plot(cvaldvdc,keffE02(31,:),'linewidth',2)
plot(cvaldvdc,keff0E02(31,:),'--')
%plot(cvaldvdc,keff(31,:),'linewidth',2)
%plot(cvaldvdc,keff2(31,:),'linewidth',2)
%plot(cvaldvdc,keff0(31,:),'--')
axis([0.25 0.75 0 0.25])
legend('E=0.02,k=10^{-6},x_0=-0.02','E=0.02,k=10^{-6},x_0=0','E=0.02,k=10^{-4},x_0=-0.02','E=0.02,k=10^{-4},x_0=0','E=0.125,k=10^{-4},x_0=-0.02','E=0.125,k=10^{-4},x_0=-0.16','E=0.125,k=10^{-4},x_0=0')
set(gca,'fontsize',22)
xlabel('dye concentration','fontsize',22)
%add arrows by hand using insert arrow
ylabel('\kappa_{eff}','fontsize',22)
%% E0125 cross-sections using cvaldvdc
%load('E0125compare3v2.mat')
figure; contourf(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(cAllnorm(50,:,:,5)),cvaldvdc); shading 'flat'
caxis([0 1]); xlabel('x','fontsize',14); ylabel('z','fontsize',14); title('C, E=0.125, x_0=-0.02, k=10^{-4}, t=39','fontsize',14); colorbar
set(gca,'fontsize',12)

figure; contourf(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(cAllnorm0(50,:,:,5)),cvaldvdc)
caxis([0 1]); xlabel('x','fontsize',14); ylabel('z','fontsize',14); title('C, E=0.125, x_0=0, k=10^{-4}, t=39','fontsize',14); colorbar
set(gca,'fontsize',12)

figure; contourf(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(cAllnorm2(50,:,:,5)),cvaldvdc)
caxis([0 1]); xlabel('x','fontsize',14); ylabel('z','fontsize',14); title('C, E=0.125, x_0=-0.16, k=10^{-4}, t=39','fontsize',14); colorbar
set(gca,'fontsize',12)

figure; contourf(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(cAllnorm(50,:,:,31)),cvaldvdc)%0.46:0.005:0.56%,0.1:0.02:0.57)
caxis([0.46 0.56]); 
xlabel('x','fontsize',14); ylabel('z','fontsize',14); title('C, E=0.125, x_0=-0.02, k=10^{-4}, t=299','fontsize',14); colorbar
set(gca,'fontsize',12)

figure; contourf(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(cAllnorm0(50,:,:,31)),cvaldvdc)%0.46:0.005:0.56%,0.15:0.02:0.57)
caxis([0.46 0.56]); 
xlabel('x','fontsize',14); ylabel('z','fontsize',14); title('C, E=0.125, x_0=0, k=10^{-4}, t=299','fontsize',14); colorbar
set(gca,'fontsize',12)

figure; contourf(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(cAllnorm2(50,:,:,31)),cvaldvdc)%0.46:0.005:0.56%,0.15:0.02:0.57)
caxis([0.46 0.56]); 
xlabel('x','fontsize',14); ylabel('z','fontsize',14); title('C, E=0.125, x_0=-0.16, k=10^{-4}, t=299','fontsize',14); colorbar
set(gca,'fontsize',12)
%%
figure; contourf(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(keffmap(50,:,:,5)),sort(unique(keff(5,:))))%,0:0.002:0.016)%,0:0.05:0.27)
caxis([0 0.015]); shading 'flat'
xlabel('x','fontsize',14); ylabel('z','fontsize',14); title('K, E=0.125, x_0=-0.02, k=10^{-4}, t=39','fontsize',14); colorbar
set(gca,'fontsize',12)

figure; pcolor(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(keffmap0(50,:,:,5)))%,0:0.003:0.016)%,0:0.05:0.27)
caxis([0 0.015]); shading 'flat'
xlabel('x','fontsize',14); ylabel('z','fontsize',14); title('K, E=0.125, x_0=0, k=10^{-4}, t=39','fontsize',14); colorbar
set(gca,'fontsize',12)

figure; pcolor(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(keffmap2(50,:,:,5)))%,0:0.003:0.016)
caxis([0 0.015]); shading 'flat'
xlabel('x','fontsize',14); ylabel('z','fontsize',14); title('K, E=0.125, x_0=-0.16, k=10^{-4}, t=39','fontsize',14); colorbar
set(gca,'fontsize',12)

figure; pcolor(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(keffmap(50,:,:,31)))%,[0 0.0001 0.00025 0.001:0.002:0.014])%,[0 0.005 0.01:0.01:0.1])
caxis([0 0.016]); shading 'flat'
xlabel('x','fontsize',14); ylabel('z','fontsize',14); title('K, E=0.125, x_0=-0.02, k=10^{-4}, t=299','fontsize',14); colorbar
set(gca,'fontsize',12)

figure; pcolor(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(keffmap0(50,:,:,31)))%,[0 0.0001 0.00025 0.001:0.002:0.014])%,[0 0.005 0.01:0.01:0.1])
caxis([0 0.016]); shading 'flat'
xlabel('x','fontsize',14); ylabel('z','fontsize',14); title('K, E=0.125, x_0=0, k=10^{-4}, t=299','fontsize',14); colorbar
set(gca,'fontsize',12)

figure; pcolor(squeeze(x(50,:,:)),squeeze(z(50,:,:)),squeeze(keffmap2(50,:,:,31)))%,[0 0.0001 0.00025 0.001:0.002:0.014])%%caxis([0 0.025]); 
caxis([0 0.016]); shading 'flat'
xlabel('x','fontsize',14); ylabel('z','fontsize',14); title('K, E=0.125, x_0=-0.16, k=10^{-4}, t=299','fontsize',14); colorbar
set(gca,'fontsize',12)

figure; plot(time,keffInt,'LineWidth',2); hold all; plot(time,keffInt0,'--','LineWidth',2); plot(time,keffInt2,'LineWidth',2); 
plot(time,ones(size(time))*(pi^6)*D_val./8,'k--'); plot(time,ones(size(time))*8*(pi^3)*D_val,'r--')
legend('x0=-0.02','x0=0','x0=-0.16','t\rightarrow\infty, ct','t\rightarrow\infty, st');
xlabel('time','fontsize',14); ylabel('\int K_{eff} dV','fontsize',14); title('E=0.125, k=10^{-4}','fontsize',14)
set(gca,'fontsize',12); axis tight
 