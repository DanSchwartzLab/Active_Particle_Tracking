%% Inverse opal project

% Written by: Haichao Wu
% Adapted by: Anni Shi
% Date: 8/22/2017
% Updated date: 11/3/2022

% Algorithm:
% 1. Find the maximum and minimum of a stage
% 2. If the absolute value of the difference of next position value and 
% maximum/minimum is less than a certain threshold (e.g. 500), assume this
% is a jump, then set the initial value as max and min, and update again
% 3. Need to consider x, y and z together, as long as one of them change
% into next stage, update max and min for all of these three.


%% Clear data


clc
clear
%% Part 1: Find the jump  frequency and dwell time
%Import the data
filename=dir('*.csv');
filename={filename.name};


DwellT=[]; %To save the dwelling time information

for jj=1:length(filename)
    data=xlsread(filename{jj});
    t=data(:,1);
    P=data(:,2:4);
%stage_change consider x,y,z together
x=P(:,1); y=P(:,2);z=P(:,3);
% for x part
max_x=x(1);
min_x=x(1);
max_y=y(1);
min_y=y(1);
max_z=z(1);
min_z=z(1);
n=numel(t);
m=1; % stage number
T=t(1);
for i=2:n
    if abs(max_x-x(i))<500&&abs(min_x-x(i))<500
        if x(i)>=max_x
            max_x=x(i);
        end
        if x(i)<=min_x
            min_x=x(i);
        end
    else
        tt(m)=t(i);
        ind(m)=i;
        m=m+1;
        T=t(i);
        max_x=x(i);
        min_x=x(i);
        
        % also need to update the y & z value
        max_y=y(i);  
        min_y=y(i);
        
        max_z=z(i);
        min_z=z(i);
    end
    
    if abs(max_y-y(i))<500&&abs(min_y-y(i))<500
        if y(i)>=max_y
            max_y=y(i);
        end
        if y(i)<=min_y
            min_y=y(i);
        end
    else
        tt(m)=t(i);
        ind(m)=i;
        m=m+1;
        T=t(i);
        max_y=y(i);
        min_y=y(i);
        
        max_x=x(i);
        min_x=x(i);
        max_z=z(i);
        min_z=z(i);
    end
   
    if abs(max_z-z(i))<515&&abs(min_z-z(i))<515  % z has larger threshold due to the drifting effect
        if z(i)>=max_z
            max_z=z(i);
        end
        if z(i)<=min_z
            min_z=z(i);
        end
    else
        tt(m)=t(i);
        ind(m)=i;
        m=m+1;
        T=t(i);
        max_z=z(i);
        min_z=z(i);
        
        max_x=x(i);
        min_x=x(i);
        max_y=y(i);
        min_y=y(i);
        
    end
    
end




ACtime=0.05;   % acquisition time 50ms
for i=1:length(ind)-1
    deltaT(i)=(ind(i+1)-ind(i))*0.05;
end

JumpF=m-1;
F(jj)=JumpF;%number of jumps for each file
DwellT=cat(1,DwellT,deltaT');

figure();
plot(t,x);hold on;plot(t,y);plot(t,z);
for i=1:length(ind)
    xd=[t(ind),t(ind)];yd=[-2000,9000];
    plot(xd,yd,'k');
end
xlabel('Time(s)')
ylabel('Position(nm)')


ind=[]; %clear the data
deltaT=[];

end

%% Find the Dwell time CDF
E=[0:1:50];
[nn,XX]=histcounts(DwellT,E,'Normalization','cdf');
size=size(DwellT);
Num=size(1);
xx=[0,XX(1:end-1)+diff(XX)/2];
nn=[0,nn];
mm=1-nn;
std=sqrt(mm*Num*1)/Num;
figure()
plot(xx,mm,'o-');
xlabel('dwell time(s)');
ylabel('frequency')
DATA=[xx',mm',std'];

%% Probability of wall-particle distance
S=[0:5:180];
[dp,S] = histcounts(wallD,S,'Normalization','probability');
ss=S(2:end);
distData=[ss',dp'];
scatter(ss',dp')
%% MSD Analyzer
%Need @msdnalyzer in path
clc
clear

filename=dir('*.csv');
filename={filename.name};
NumberTrajectories = length(filename);

for i=1:1:length(filename);
    data=xlsread(filename{i});
      time=data(:,1)*0.01;% * step time: 10 ms
      X=data(:,2:3)/1000;%convert to µm
%     time=data(:,1);
%     X=data(:,2:3);  
    tracks{i} = [time X];
    
end

% Initialize MSDanalyzer
ma = msdanalyzer(2, 'um', 's');
ma = ma.addAll(tracks);
mmsd = ma.getMeanMSD;
%% mean msd
mmsd = ma.getMeanMSD;
t = mmsd(:,1);
x = mmsd(:,2);
dx = mmsd(:,3) ./ sqrt(mmsd(:,4));
errorbar(t, x, dx, 'k');
%% 
% Plot MSd curves
figure
ma.plotMSD;
legend;


