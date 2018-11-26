%using boxplot to put together comparison of ranges of FTLES and FALS
%for oral exam presentation
%taken from http://www.ics.uci.edu/~vpalepu/2013/03/03/Plotting-Box-plots-in-Groups-for-Vectors-of-Varying-Lengths.html

%Resonant Layers
%resonant layer width range on left box ('c'), resonant layer FALS range on right
%box ('y') in each pair
load('fals2.mat')
U=0.1;
% lamdaMinR=ftleResMin.*U./D;
% lamdaMaxR=ftleResMax.*U./D;
% lamdaMinC=ftleChaosMainMin.*U./D;
% lamdaMaxC=ftleChaosMainMax.*U./D;

lamdaMinR=ftleResMin2.*U./D;
lamdaMaxR=ftleResMax2.*U./D;
lamdaMinC=ftleChaosMainMin2.*U./D;
lamdaMaxC=ftleChaosMainMax2.*U./D;

falsDimRmax=(0.0103./(lamdaMinR)).^(1/0.85);%in cm
falsNondimRmax=falsDimRmax./(100.*D);
falsDimRmin=(0.0103./(lamdaMaxR)).^(1/0.85);%in cm
falsNondimRmin=falsDimRmin./(100.*D);
falsDimCmax=(0.0103./(lamdaMinC)).^(1/0.85);%in cm
falsNondimCmax=falsDimCmax./(100.*D);
falsDimCmin=(0.0103./(lamdaMaxC)).^(1/0.85);%in cm
falsNondimCmin=falsDimCmin./(100.*D);

resWidthMax(1)=0.04;
resBoxes=[falsNondimRmin;falsNondimRmax;resWidthMin;resWidthMax];
chaosBoxes=[falsNondimCmin;falsNondimCmax;chaosMainWidthMin;chaosMainWidthMax];
%%
x = resBoxes(:,1:4);
y = chaosBoxes(:,1:4); %[0.01,0.04,0.014,0.066,0.01,0.03,0.012,0.017,0.02,0.07,0.015,0.023,0.01,0.04,0.019,0.027];
z=cat(2,x,y);
group = [1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13,14,14,15,15,16,16];
positions = [1 1.5 2.5 3 4.5 5 6 6.5 7.5 8 9 9.5 10.5 11 12 12.5];
figure
boxplot(z,group, 'positions', positions);
color = ['c','y', 'c','y', 'c','y', 'c', 'y','c','y', 'c','y', 'c','y', 'c', 'y'];
h = findobj(gca,'Tag','Box');
 for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),color(j),'FaceAlpha',.5);
 end
 title('Layer Width (blue) vs FALS (yellow)','fontsize',16)
 set(gca,'XTick',[1.25 2.75 4.75 6.25 7.75 9.25 10.75 12.25]);
 set(gca,'XTickLabel',{'E 0.25', 'E 0.125','E 0.02','E 0.0005','E 0.25', 'E 0.125','E 0.02','E 0.0005'});
 xlabel('Left: resonant region. Right:outer region.','fontsize',14)
 ylabel('Widths','fontsize',14)
 ylim([0 0.3])
 %%
 %Main chaotic region
%main chaotic region width range on left box ('c'), FALS range on right
%box ('y') in each pair
% x = chaosBoxes(:);%[0.25,0.35,0.01,0.017,0.25,0.35,0.01,0.022,0.1,0.2,0.012,0.023,0.05,0.15,0.013,0.019];
% group = [1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8];
% positions = [1 1.5 2.5 3 4 4.5 5.5 6];
% figure(2)
% boxplot(x,group, 'positions', positions);
% color = [ 'c','y', 'c','y', 'c', 'y', 'c','y'];
% h = findobj(gca,'Tag','Box');
%  for j=1:length(h)
%     patch(get(h(j),'XData'),get(h(j),'YData'),color(j),'FaceAlpha',.5);
%  end
%  title('Outer Chaotic Region Width (blue) vs FALS (yellow)','fontsize',16)
%   xlabel('Left to right E=0.25 0.125 0.02 0.0005','fontsize',14)
%    ylabel('Widths','fontsize',14)
   
kappaRmin=0.0103.*(falsDimRmin).^1.15;
kappaRmax=0.0103.*(falsDimRmax).^1.15;
kappaCmin=0.0103.*(falsDimCmin).^1.15;
kappaCmax=0.0103.*(falsDimCmax).^1.15;
tauRmin=((D.*resWidthMin).^2)./(kappaRmax./(1e4));
tauRmax=((D.*resWidthMin).^2)./(kappaRmin./(1e4));
tauCmin=((D.*resWidthMin).^2)./(kappaCmax./(1e4));
tauCmax=((D.*resWidthMin).^2)./(kappaCmin./(1e4));

tauRmin2=((100*D.*resWidthMin).^0.85)./(3600*0.0103);%uses tau=L^2/kappa and okubo kappa(cm^2/s)=0.0103L(cm)^1.15 
tauRmax2=((100*D.*resWidthMax).^0.85)./(3600*0.0103);%tau in hours!
kappa2min=0.0103*((100*D.*resWidthMin).^1.15)/1e4;
kappa2max=0.0103*((100*D.*resWidthMax).^1.15)/1e4;%okubo, 1e4 converts back to m^2/s
%tauCmin2=((D.*resWidthMin).^2)./(kappaCmax./(1e4));
%tauCmax2=((D.*resWidthMin).^2)./(kappaCmin./(1e4));

tauBoxes=[tauRmin2;tauRmax2];%tauCmin;tauCmax];
x=tauBoxes(:,2:3);%/3600;
%group = [1 1 2 2];
%positions = [1 1.5 2.5 3];
figure
boxplot(x);%,group, 'positions', positions);
color = ['w','w','w','w'];%'b','w','b','w'];
h = findobj(gca,'Tag','Box');
 for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),color(j),'FaceAlpha',.5);
 end
 title('Diffusive Island Crossing Time','fontsize',16)
 %xlabel('Left to right E=0.125 0.02','fontsize',14)
 ylabel('Time in Hours','fontsize',14)
%  set(gca,'XTick',[1 2 3 4]);
%  set(gca,'XTickLabel',{'E 0.25','E 0.125','E 0.02','E 0.0005'});
%  ylim([0 370])
  set(gca,'XTick',[1 2]);
 set(gca,'XTickLabel',{'E 0.125','E 0.02'});
 ylim([0 10])
%%
%% NEW Version
load('fals3.mat')

%% nondimensional filament arrest width and region width
resBoxes=[deltaNDmin;deltaNDmax;resWidthMin;resWidthMax];
chaosBoxes=[deltaNDmin;deltaNDmax;chaosMainWidthMin;chaosMainWidthMax];
x = resBoxes(:,1:4);
y = chaosBoxes(:,1:4); %[0.01,0.04,0.014,0.066,0.01,0.03,0.012,0.017,0.02,0.07,0.015,0.023,0.01,0.04,0.019,0.027];
z=cat(2,x,y);
group = [1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13,14,14,15,15,16,16];
positions = [1 1.5 2.5 3 4.5 5 6 6.5 7.5 8 9 9.5 10.5 11 12 12.5];
figure
boxplot(z,group, 'positions', positions);
color = ['c','y', 'c','y', 'c','y', 'c', 'y','c','y', 'c','y', 'c','y', 'c', 'y'];
h = findobj(gca,'Tag','Box');
 for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),color(j),'FaceAlpha',.5);
 end
 title('Layer Width (blue) vs \delta (yellow)','fontsize',16)
 set(gca,'XTick',[1.25 2.75 4.75 6.25 7.75 9.25 10.75 12.25]);
 set(gca,'XTickLabel',{'E 0.25', 'E 0.125','E 0.02','E 0.0005','E 0.25', 'E 0.125','E 0.02','E 0.0005'});
 xlabel('Left: resonant region. Right:outer region.','fontsize',14)
 ylabel('Nondimensional Widths','fontsize',14)
 ylim([0 0.3])
 hold on; plot([7 7],[0 0.3],'k')
 %% dimensional filament arrest width and region width
resBoxes=[deltaCMmin;deltaCMmax;resCMmin;resCMmax]./100;
chaosBoxes=[deltaCMmin;deltaCMmax;chaosMainCMmin;chaosMainCMmax]./100;
x = resBoxes(:,1:4);
y = chaosBoxes(:,1:4); %[0.01,0.04,0.014,0.066,0.01,0.03,0.012,0.017,0.02,0.07,0.015,0.023,0.01,0.04,0.019,0.027];
z=cat(2,x,y);
group = [1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13,14,14,15,15,16,16];
positions = [1 1.5 2.5 3 4.5 5 6 6.5 7.5 8 9 9.5 10.5 11 12 12.5];
figure
boxplot(z,group, 'positions', positions);
color = ['c','y', 'c','y', 'c','y', 'c', 'y','c','y', 'c','y', 'c','y', 'c', 'y'];
h = findobj(gca,'Tag','Box');
 for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),color(j),'FaceAlpha',.5);
 end
 title('Layer Width (blue) vs \delta (yellow)','fontsize',16)
 set(gca,'XTick',[1.25 2.75 4.75 6.25 7.75 9.25 10.75 12.25]);
 set(gca,'XTickLabel',{'E 0.25', 'E 0.125','E 0.02','E 0.0005','E 0.25', 'E 0.125','E 0.02','E 0.0005'});
 xlabel('Left: resonant region. Right:outer region.','fontsize',14)
 ylabel('Widths in m','fontsize',14)
 set(gca,'YScale','log')
ylim([1e-1 1e3])
 hold on; plot([7 7],[1e-1 1e3],'k')
 %% diffusive island crossing time vs contracting ftle time
 
resBoxes=[1./lambda3maxDim;1./lambda3minDim;tauResMin;tauResMax]./3600;
%chaosBoxes=[deltaCMmin;deltaCMmax;chaosMainCMmin;chaosMainCMmax];
x = resBoxes(:,1:4);
%y = chaosBoxes(:,1:4); %[0.01,0.04,0.014,0.066,0.01,0.03,0.012,0.017,0.02,0.07,0.015,0.023,0.01,0.04,0.019,0.027];
%z=cat(2,x,y);
group = [1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8];%,9,9,10,10,11,11,12,12,13,13,14,14,15,15,16,16];
positions = [1 1.5 2.5 3 4.5 5 6 6.5];% 7.5 8 9 9.5 10.5 11 12 12.5];
figure
boxplot(x,group, 'positions', positions);
color = ['c','y', 'c','y', 'c','y', 'c', 'y'];%,'c','y', 'c','y', 'c','y', 'c', 'y'];
h = findobj(gca,'Tag','Box');
 for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),color(j),'FaceAlpha',.5);
 end
 title('Diffusive (blue) vs Advective (yellow) Timescale','fontsize',16)
 set(gca,'XTick',[1.25 2.75 4.75 6.25]);% 7.75 9.25 10.75 12.25]);
 set(gca,'XTickLabel',{'E 0.25', 'E 0.125','E 0.02','E 0.0005'});%,'E 0.25', 'E 0.125','E 0.02','E 0.0005'});
 %xlabel('Left: resonant region. Right:outer region.','fontsize',14)
 ylabel('Times in hours','fontsize',14)
 ylim([0 100])
 
 %%