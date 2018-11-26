function dqdt = TurbulentRotatingCanEq2(t,q,dt,R,E,eps,pert,x0,sigma)
%recently corrected the polar to cartesian conversion!
%
%t=time
%q=initial position(s)
%dt=spacing of t, assumed constant
%R=radius of the can, nondim a in the manuscript
%eps=strength of velocity perturbations
%pert is the perturbation form option. 1 has no w perturbation, 2 is the 
%2-function 3d perturbation. 4 is Irina's 2d perturbation.
%x0=offset for perturbation
%sigma=turbulence parameter, a constant scalar
%sigma=0.01 gives mixing ratio of 10^-13 through at least t=10k


%%%Constants
n=length(q);
timeIndex=t/dt+1;%used with turbulence
Eratio=sqrt(E);%ratio of Ekman layer thickness to dimensional tank depth
s=sin(0.5/Eratio);
c=cos(0.5/Eratio);
s2=sinh(0.5/Eratio);
c2=cosh(0.5/Eratio);
A=-0.5*c*s2/((s^2)*(c2^2)+(c^2)*(s2^2));
B=0.5*c2*s/((s^2)*(c2^2)+(c^2)*(s2^2));
D=A*(s*s2-c*c2)+B*(c*c2+s*s2);

%%%%Variables
xx=q(1:3:end-2); 
yy=q(2:3:end-1);  
zz=q(3:3:end);

[theta, r]=cart2pol(xx,yy);
r2=sqrt((yy.^2)+(xx-x0).^2);
zeta=(zz-0.5)/Eratio;%modified vertical variable

%%%%Calculate velocities
ur=r.*((R-r).^2).*(A.*sin(zeta).*cosh(zeta)+B.*cos(zeta).*sinh(zeta));
vtheta=r.*((R-r).^2).*(0.5+B.*sin(zeta).*cosh(zeta)-A.*cos(zeta).*sinh(zeta));
%vtheta could be multiplied by a constant without affecting the continuity or boundary conditions
w=-Eratio.*(R-r).*(R-2.*r)...
        .*(A.*(sin(zeta).*sinh(zeta)-cos(zeta).*cosh(zeta))...
        +B.*(sin(zeta).*sinh(zeta)+cos(zeta).*cosh(zeta))-D);

u=cos(theta).*ur-r.*sin(theta).*vtheta;
v=sin(theta).*ur+r.*cos(theta).*vtheta;
    
switch pert
    case 1
        u=u+eps.*sinh(zz./Eratio).*(R-r).*((-yy.*(r2.^2)./r)+(yy.*(R-r)))./sinh(1/Eratio);
        v=v+eps.*sinh(zz./Eratio).*(R-r).*((-(xx-x0).*(R-r))+(xx.*(r2.^2)./r))./sinh(1/Eratio);
    case 2
        u=u+eps.*(-0.25).*r2.*yy.*((R-r).^2).*(pi.*zz.*cos(pi.*zz)+sin(pi.*zz));
        v=v+eps.*0.25.*r2.*(xx-x0).*((R-r).^2).*(pi.*zz.*cos(pi.*zz)+sin(pi.*zz));
        w=w+eps.*0.5.*r2.*x0.*yy.*(R-r).*zz.*sin(pi.*zz)./r;
    case 4
        u=5.*u+eps.*(2.*yy.*(R - xx.^2 - yy.^2) + 2.*yy.*(4*R^2 - (-x0 + xx).^2 -... 
            yy.^2)).*(sinh(zz./Eratio)./sinh(1./Eratio));
        v=5.*v+eps.*(-2.*(-x0 + xx).*(R - xx.^2 - yy.^2) - 2.*xx.*(4*R^2 - (-x0 + ...
            xx).^2 - yy.^2)).*(sinh(zz./Eratio)./sinh(1./Eratio));
        w=5.*w;
    case 5 %2d, modified Irina's to get correct BCs
        u=u+eps.*(R - r.^2).*yy.*(2.*(R-r.^2).* + 4.*(4*R^2 - r2.^2))... 
        .*(sinh(zz./Eratio)./sinh(1./Eratio));
        v=v+eps.*(R - r.^2).*(-2.*(-x0 + xx).*(R - r.^2) - 4.*xx.*(4*R^2 - r2.^2)) ...
        .*(sinh(zz./Eratio)./sinh(1./Eratio));
        %w=w;
    case 6 %3d, based off 4
        u=5.*u+eps.*(R-r.^2).*(4*R^2 - r2.^2).*yy.*(0.5.*sin(pi.*zz).*(R-r.^2)-zz.*cos(pi.*zz).*(4*R^2 - r2.^2));
        v=5.*v+eps.*(R-r.^2).*(4*R^2 - r2.^2).*0.5.*(xx.*zz.*cos(pi.*zz).*(4*R^2 - r2.^2)-sin(pi.*zz).*(xx-x0).*(R-r.^2));
        w=5.*w+eps.*2.*(R-r.^2).*(4.*(R.^2)-r2.^2).*x0.*zz.*sin(pi.*zz);
end

if sigma~=0 && mod(timeIndex,1)==0 %used for turbulence cases before
    u=u+sigma*randn(size(u));
    v=v+sigma*randn(size(v));
    w=w+sigma*randn(size(w));
end


%nout=find(((xx.^2+yy.^2)>(R^2))|(zz<0)|(zz>1)); u(nout)=0; v(nout)=0; w(nout)=0; % set velocity outside the can to 0
%nouta=find(((xx.^2)>(R^2))|((yy.^2)>(R^2))|(zz<0)|(zz>1)); u(nouta)=0; v(nouta)=0; w(nouta)=0;

%nudge back inside cylinder
u(((xx.^2+yy.^2)>(R^2)) & xx>0)=-0.001;
u(((xx.^2+yy.^2)>(R^2)) & xx<0)=0.001;
v(((xx.^2+yy.^2)>(R^2)) & yy>0)=-0.001;
v(((xx.^2+yy.^2)>(R^2)) & yy<0)=0.001;
w(zz>1)=-0.001;
w(zz<0)=0.001;
%%%%return
dqdt=[u; v; w];
dqdt=reshape(dqdt,[n/3 3])';
dqdt=dqdt(:); 