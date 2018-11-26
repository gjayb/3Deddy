function [ dqdt ] = ellipsoidEq( t,q,lambda,ubar )
%ellipsoidEq linear strain velocity field
%   lambda, ubar 1x3, t not really needed, q 1x3n
n=length(q);
xx=q(1:3:end); yy=q(2:3:end); zz=q(3:3:end);
dxdt=ubar(1)+lambda(1)*xx;
dydt=ubar(2)+lambda(2)*yy;
dzdt=ubar(3)+lambda(3)*zz;

dqdt=[dxdt; dydt; dzdt];
dqdt=reshape(dqdt,[n/3 3])';
dqdt=dqdt(:); 
end

