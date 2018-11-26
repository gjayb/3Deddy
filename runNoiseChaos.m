%running ensembles of trajectories to compare stochastic and steady chaotic
%perturbations in the rotating can.
clear
R=1;
pert1=4;
x01=-0.5;

dt=0.01;%time spacing for calculation; needs to be small
t=0:dt:10000;%time
dtrecord=1;%time spacing for output; larger than dt to not have giant vectors
ndt=dtrecord/dt;
%timesWanted=0:dtrecord:max(t);
timesIndices=1:ndt:length(t);

%Eks=[0.125 0.02 0.01 0.005 0.002 0.0005];
%nzInner=[1 1 1 2 2 1];
%zInner=[NaN NaN NaN 0.8078 0.8779 NaN];
Eks=[0.25 0.125 0.02 0.0005];

load 'initialSpheres.mat'

for Eki=3%2:3
    E=Eks(Eki)
    
    for rInit=[0 0.1 0.4]%0:0.1:0.5%[0.1 0.3 0.7 0.9]%[0.5]% 0.9] 
        if rInit==0.5
             qin=[x1';y1';z1'];
        else
            dr=rInit-0.5;
            %z22=z2+dz;
            %qin=[x2';y2';z22'];
            y12=y1+dr;
            qin=[x1';y12';z1'];
        end
        for pertType=8:9%1:6%[1 3]
            switch pertType
                case 1
                    sigma=0.042
                    eps=0.0
                case 2
                    sigma=0.0
                    eps=0.01
                case 3
                    sigma=0.042
                    eps=0.01
                case 4
                    sigma=0.013
                    eps=0.0
                case 5
                    sigma=0.13
                    eps=0.0
                case 6
                    sigma=0
                    eps=0.08
                case 7
                    sigma=0.013
                    eps=0.01
                case 8
                    sigma=0.013
                    eps=0.08
                case 9
                    sigma=0.042
                    eps=0.08
            end

           
            qq=ode4(@TurbulentRotatingCanEq2,t,qin(:),dt,R,E,eps,pert1,x01,sigma);
            xtr=qq(timesIndices,1:3:end-2);
            ytr=qq(timesIndices,2:3:end-1);
            ztr=qq(timesIndices,3:3:end);
            ttr=t(timesIndices);
            clear qq;
            [theta,r]=cart2pol(xtr,ytr);

            ftitle1s=strcat('E',num2str(E),'sigma',num2str(sigma),'eps',num2str(eps),'r',num2str(rInit),'S.mat');
            save(ftitle1s)
        end %for pert
    end %for rInit
    
%     if ~isnan(zInner(Eki))
%         zInit=zInner(Eki);
%         [thetaIn,rIn,zIn]=ndgrid(thetaInit,rInit,zInit);
%         [xIn,yIn]=pol2cart(thetaIn,rIn);
%         for pertType=1:3
%             switch pertType
%                 case 1
%                     sigma=0.013;
%                     eps=0.0;
%                 case 2
%                     sigma=0.0;
%                     eps=0.01;
%                 case 3
%                     sigma=0.0065;
%                     eps=0.005;
%             end
% 
%             qin=[xIn(:).';yIn(:).';zIn(:).'];
%             qq=ode4(@TurbulentRotatingCanEq,t,qin(:),dt,R,E,eps,pert1,x01,sigma);
%             xtr=qq(timesIndices,1:3:end-2);
%             ytr=qq(timesIndices,2:3:end-1);
%             ztr=qq(timesIndices,3:3:end);
%             ttr=t(timesIndices);
%             clear qq;
%             clear t;
%             [theta,r]=cart2pol(xtr,ytr);
% 
%             ftitle1s=strcat('E',num2str(E),'sigma',num2str(sigma),'eps',num2str(eps),'r',num2str(rInit),'z',num2str(zInit),'.mat');
%             save(ftitle1s)
%         end
%     end
end %for E
    