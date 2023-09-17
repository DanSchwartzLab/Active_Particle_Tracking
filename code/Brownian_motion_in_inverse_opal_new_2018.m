% Brownian motion in a cavity
% Written by: Haichao Wu
% Date: 11/8/2017


%{
Important info from experiments:
According to 45nm particle, 120nm channel size video, for 50ms, mean displacement is 
0nm, SD is 67nm. For the z direction, could be a little smaller, this is 
due to experimental error induced by tracking method.
%}

%{
steps for run this simulation
1. Change channel size, particle size, reps, max_steps
2. Change the name of three figures (MSD, CDF, fitted_CDF)
3. Change the name of saved workspace
After finish running:
1. Make sure figure and workspace is appropriately saved, and copy these file
into a specific folder
2. Copy the fitting results in the command window
%}
clc
clear

tic

reps=500; %Number of trajectory simulated
Cr=250; %Cavity radius
Pr=46/2; %Tracer particle radius  
R=Cr-Pr;  %Accesible radius
r=70/2;  %radius of connecting channel
q=asin((r-Pr)/R); %half of the effective opening angle 
max_steps=20000;
N=zeros(max_steps,3);
N_step=zeros(reps,1); %save the number of steps a particle need to jump out of cavity

Num=0; % How many times the particle crash on the wall.

sigma0=95;
mu=0; sigma=80; %Get from experimental data for generating step size

%data space for save MSD info
rx2=zeros(reps,max_steps); %vector to store squared displacements
stdevx=zeros(reps,max_steps); %vector to store error in squared displacement
numx=zeros(reps,max_steps); %number of steps for each interval

for h=1:reps
    
% update the data for new trajectory
x=zeros(max_steps,1);
y=zeros(max_steps,1);
z=zeros(max_steps,1);
sum1=zeros(max_steps,1);
x=single(x);
y=single(y);
z=single(z);

%set up the initial position of the tracer particle
x(1)=normrnd(mu,sigma0);
y(1)=normrnd(mu,sigma0);
z(1)=normrnd(mu,sigma0);

    
for i=2:max_steps
    x(i)=x(i-1)+normrnd(mu,sigma);
    y(i)=y(i-1)+normrnd(mu,sigma);
    z(i)=z(i-1)+normrnd(mu,sigma);
    N(i,:)=[x(i),y(i),z(i)];
    sum1(i)=sqrt(x(i)^2+y(i)^2+z(i)^2);
    if norm(N(i,:))>R  % test if the point is beyond the cavity
        %Transform Cartesian coordinates to spherical(Azimuth and Elevation Angles)
        [azi,ele,rad]=cart2sph(x(i),y(i),z(i)); 
        % constraints on the horizontal direction
        if (azi>=-q)&&(azi<=q)&&(ele>=-q)&&(ele<=q)
            N_step(h)=i;
            disp('success jump');
            break;
        elseif (azi>=pi/3-q)&&(azi<=pi/3+q)&&(ele>=-q)&&(ele<=q)
            N_step(h)=i;
            disp('success jump');
            break;        
        elseif (azi>=2*pi/3-q)&&(azi<=2*pi/3+q)&&(ele>=-q)&&(ele<=q)
            N_step(h)=i;
            disp('success jump');
            break;
        elseif (azi>=pi-q)&&(azi<=pi)&&(ele>=-q)&&(ele<=q)
            N_step(h)=i;
            disp('success jump');
            break;
        elseif (azi>=-pi)&&(azi<=-pi+q)&&(ele>=-q)&&(ele<=q)
            N_step(h)=i;
            disp('success jump');
            break;
        elseif (azi>=-2*pi/3-q)&&(azi<=-2*pi/3+q)&&(ele>=-q)&&(ele<=q)
            N_step(h)=i;
            disp('success jump');
            break;
        elseif (azi>=-pi/3-q)&&(azi<=-pi/3+q)&&(ele>=-q)&&(ele<=q)
            N_step(h)=i;
            disp('success jump');
            break;
        % constraints not on the horizontal direction
        elseif (azi>=pi/6-q)&&(azi<=pi/6+q)&&(ele>=pi/4-q)&&(ele<=pi/4+q)
            N_step(h)=i;
            disp('success jump');
            break;
        elseif (azi>=5*pi/6-q)&&(azi<=5*pi/6+q)&&(ele>=pi/4-q)&&(ele<=pi/4+q)
            N_step(h)=i;
            disp('success jump');
            break;
        elseif (azi>=-pi/2-q)&&(azi<=-pi/2+q)&&(ele>=pi/4-q)&&(ele<=pi/4+q)
            N_step(h)=i;
            disp('success jump');
            break;
        elseif (azi>=pi/2-q)&&(azi<=pi/2+q)&&(ele>=-pi/4-q)&&(ele<=-pi/4+q)
            N_step(h)=i;
            disp('success jump');
            break;
        elseif (azi>=-5*pi/6-q)&&(azi<=-5*pi/6+q)&(ele>=-pi/4-q)&&(ele<=-pi/4+q)
            N_step(h)=i;
            disp('success jump');
            break;
        elseif (azi>=-pi/6-q)&&(azi<=-pi/6+q)&(ele>=-pi/4-q)&&(ele<=-pi/4+q)
            N_step(h)=i;
            disp('success jump');
            break;
        end
        
        a0=N(i-1,:);
        a1=N(i,:);

        d=norm(a1-a0);
        A0=[a1+1,1];
        f=@(x) myfun(x,a0,a1,R);
        A=fsolve(f,A0);
        A1=[A(1),A(2),A(3)]; %calculate the intersection point
        
        

        d1=norm(A1-a0); % calculate the  Euclidean distance
        D=d-d1;
        f2=@(x) myfun2(x,a0,A1,D);
        B0=[0,0,0,1];
        lb=[-R,-R,-R,-1000]; %constraints lower bound
        ub=[R,R,R,1000]; %constraints upper bound
        con=@(x) myfun3(x,R);
        B=fmincon(f2,B0,[],[],[],[],lb,ub,con);
        B1=[B(1),B(2),B(3)];
        D2=norm(B1-A1);
        x(i)=B1(1);y(i)=B1(2);z(i)=B1(3);
        N(i,:)=[x(i),y(i),z(i)];
        Num=Num+1;
    end
end
% figure()


% initial radius
r0=[x(1) y(1) z(1)];
plot3(r0(1),r0(2),r0(3),'or','MarkerFaceColor','r');
hold on;

if N_step(h)~= 0
    % case 1: successfully jump
% final radius
    rf=[x(N_step(h)) y(N_step(h)) z(N_step(h))];
    plot3(x(1:N_step(h)),y(1:N_step(h)),z(1:N_step(h)))
else
    % case 2: cann't hop out of cavity
    N_step(h)=max_steps   % for the usage of calculating MSD
    rf=[x(end) y(end) z(end)]
    plot3(x,y,z)
end
plot3(rf(1),rf(2),rf(3),'ok','MarkerFaceColor','k');


% plot the cavity
center=[0,0,0];
normal1=[0,0,1];
normal2=[0,1,0];
normal3=[1,0,0];
radius=R;
plotCircle3D(center,normal1,radius);
plotCircle3D(center,normal2,radius);
plotCircle3D(center,normal3,radius);

% plot the channel in the horizontal direction
radius1=r;
center4=[242.7,0,0];
normal4=[1,0,0];
center5=[121.35,210.18,0];
normal5=[1,sqrt(3),0];
center6=[-121.35,210.18,0];
normal6=[-1,sqrt(3),0];
center7=[-242.7,0,0];
normal7=[-1,0,0];
center8=[-121.35,-210.18,0];
normal8=[-1,-sqrt(3),0];
center9=[121.35,-210.18,0];
normal9=[1,-sqrt(3),0];
plotCircle3D1(center4,normal4,radius1);
plotCircle3D1(center5,normal5,radius1);
plotCircle3D1(center6,normal6,radius1);
plotCircle3D1(center7,normal7,radius1);
plotCircle3D1(center8,normal8,radius1);
plotCircle3D1(center9,normal9,radius1);

% plot the channel in the rest direction
center10=[148.6,85.6,171.6];
normal10=[sqrt(3),1,2];
center11=[-148.6,85.6,171.6];
normal11=[-sqrt(3),1,2];
center12=[0,-171.6,171.6];
normal12=[0,-1,1];
center13=[-148.6,-85.6,-171.6];
normal13=[-sqrt(3),-1,-2];
center14=[148.6,-85.6,-171.6];
normal14=[sqrt(3),-1,-2];
center15=[0,171.6,-171.6];
normal15=[0,1,-1];
plotCircle3D1(center10,normal10,radius1);
plotCircle3D1(center11,normal11,radius1);
plotCircle3D1(center12,normal12,radius1);
plotCircle3D1(center13,normal13,radius1);
plotCircle3D1(center14,normal14,radius1);
plotCircle3D1(center15,normal15,radius1);



% Line of Sight between initial and final state
xx=linspace(r0(1),rf(1));
yy=linspace(r0(2),rf(2));
zz=linspace(r0(3),rf(3));
plot3(xx,yy,zz,'k--','LineWidth',2);
xlabel('x');ylabel('y');zlabel('z');
title('Confined diffusion in inverse opal with Hopping')
grid on;
hold off;


% Calculate MSD

for j=1:N_step(h)  %go through trajectory, calculating squared displacements
    for k=1:(N_step(h)-j)
        rx_squared=(x(k+j)-x(k))^2;
        rx2(h,j)=rx2(h,j)+rx_squared; %accumulate the sum of the SD for nn steps from all trajectories
        stdevx(h,j)=stdevx(h,j)+rx_squared*rx_squared; %accumulating the second moment of the SD for nn steps 
        numx(h,j)=numx(h,j)+1;  %keep track of how many trajectories made it this many steps; necessary for averaging later
    end
end


end

% plot MSD
Srx2=sum(rx2);
Snumx=sum(numx);
Sstdevx=sum(stdevx);
Srx2=Srx2./Snumx; %average the MSD vector
Sstdevx=Sstdevx./Snumx; %average the second moment of the SD
Sstdevx=Sstdevx-(Srx2.*Srx2); %calculate the variance of the SD
Sstdevx=sqrt(Sstdevx./Snumx);
figure();errorbar(Srx2(1:3000),Sstdevx(1:3000));
xlabel('number of steps');ylabel('MSD');
title('MSD Confined diffusion in inverse opal 45nm beads 70nm channel new 95_80nm per step_NEW1')
savefig('F1_MSD_45beads_70_new_95_80nm_per_step_NEW1')

MSDDATA=[Srx2(1:3600)',Sstdevx(1:3600)'];

% Find the Dwell time CDF
N_stepnew=N_step(find(N_step~=max_steps));
figure()
histogram(N_stepnew);
E=[0:75:20000];
[nn,XX]=histcounts(N_stepnew,E,'Normalization','cdf');
xx=[0,XX(1:end-1)+diff(XX)/2];
nn=[0,nn];
mm=1-nn;
tt=0.05*xx;
figure()
plot(tt,mm,'o-');
xlabel('residence time');
ylabel('frequency')
savefig('F2_CDF_43beads_70_new_95_80nm_NEW1.')

DATA=[xx',tt',mm'];
no
toc




%% Weighted nonlinear regression

modelFun=@(b,x) b(1).*exp(-x./b(2))+b(3).*exp(-x./b(4));
start=[0.5;2000;0.5;2000];

modelFun2=@(b,x) b(1).*exp(-(x./b(2)).^b(3));
start2=[1;2000;1];

modelFun3=@(b,x) b(1).*exp(-x./b(2));
start3=[1;2000];

modelFun0=@(b,x) b(1).*exp(-x./b(2));
start3=[1;2000];


N=numel(N_stepnew);
zs=1.96; % 95%CI z score
CI=1./(1+zs^2/N)*(2*zs.*sqrt(mm.*(1-mm)/N+zs^2/4/N^2));
w=1./CI;
w=w'

figure()
semilogy(xx,mm,'o');
ylim([0.01 1]);
xlabel('Number of steps');
ylabel('Probability')
title('Complimentary CDF 45nm beads 102nm channel new 104_40nm per step_NEWNEW2')


fitx=linspace(0,20000)';

nlm1=fitnlm(xx,mm,modelFun,start,'weight',w)
line(fitx,predict(nlm1,fitx),'color','b')


nlm2=fitnlm(xx,mm,modelFun2,start2,'weight',w)
line(fitx,predict(nlm2,fitx),'color','r')

nlm3=fitnlm(xx,mm,modelFun3,start3,'weight',w)
line(fitx,predict(nlm3,fitx),'linestyle','--','color','k');

nlm0=fitnlm(xx,mm,modelFun0,start3)
line(fitx,predict(nlm0,fitx),'color','g')


legend('Raw data','Weighted two exponential fitting','Weighted stretched exponential fitting',...
    'Weighted one exponential fitting','Unweighted one exponential fitting');
savefig('F3_CDF_fitting_45beads_102_new_114_40nm_per_step_NEWNEW2')

r = nlm2.Residuals.Raw;
figure()
% The raw residuals appear a variance proportional to the inverse of the
% weights, therefore, we standardize this variance to make the plot easier
% to interpret
plot(xx,r.*sqrt(w),'b^');
xlabel('x'); ylabel('Residuals, yFit - y');
title('Standardized residual analysis')

% SAVE WORKSPACE into a file
save('Inverse_opal_45nm_beads_120nm_channel_new_10nm_per_step_NEWNEW2.mat')




% calculate the collision point
function F=myfun(x,a0,a1,R)
F(1)=(x(1)-a0(1))-x(4)*(a1(1)-a0(1));
F(2)=(x(2)-a0(2))-x(4)*(a1(2)-a0(2));
F(3)=(x(3)-a0(3))-x(4)*(a1(3)-a0(3));
F(4)=(x(1))^2+(x(2))^2+(x(3))^2-R^2;
end


% calculate the bounce back distance
function F2=myfun2(x,a0,A1,D)
nn=A1; %the normal of the tangent plane
v=A1-a0; % the directional vector for hitting the wall
u=(dot(v,nn)/dot(nn,nn))*nn;
w=v-u;
vv=w-u; % the reflection directional vector
F2=zeros(4,1);
F2(1)=(x(1)-A1(1))-x(4)*vv(1);
F2(2)=(x(2)-A1(2))-x(4)*vv(2);
F2(3)=(x(3)-A1(3))-x(4)*vv(3);
F2(4)=sqrt((x(1)-A1(1))^2+(x(2)-A1(2))^2+(x(3)-A1(3))^2)-D;
F2=sum(F2.^2);
end

% Additional constraint for fmincon
function [c,ceq]=myfun3(x,R)
c=x(1)^2+x(2)^2+x(3)^2-R^2;
ceq=[];
end


function plotCircle3D(center,normal,radius)

theta=0:0.01:2*pi;
v=null(normal);
points=repmat(center',1,size(theta,2))+radius*(v(:,1)*cos(theta)+v(:,2)*sin(theta));
plot3(points(1,:),points(2,:),points(3,:),'r-');

end


function plotCircle3D1(center,normal,radius)

theta=0:0.01:2*pi;
v=null(normal);
points=repmat(center',1,size(theta,2))+radius*(v(:,1)*cos(theta)+v(:,2)*sin(theta));
plot3(points(1,:),points(2,:),points(3,:),'b-');

end