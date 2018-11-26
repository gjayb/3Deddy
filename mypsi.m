function psitr = mypsi(R,E,xtr,ytr,ztr)
%This function calculates the value of the unperturbed streamfunction at
%each point of the position vectors given. psi(radius/height, ekman number,
%xvalues,yvalues,zvalues) where the first two are scalars and the last
%three are vectors/matrices of the same size.

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
[~,r]=cart2pol(xtr,ytr);
zeta=(ztr-0.5)/Eratio;

%streamfunction; 
Rfun=0.5.*(r.^2).*((R-r).^2);
Ffun=(A.*(sin(zeta).*sinh(zeta)-cos(zeta).*cosh(zeta))...
    +B.*(cos(zeta).*cosh(zeta)+sin(zeta).*sinh(zeta))-D);
psitr=-Eratio.*Rfun.*Ffun;


end
