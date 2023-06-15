%Active_Particle_Code
%% 
clc
clear
close all

%% MSD Analyzer
%Need @msdnalyzer in path

filename=dir('*.csv');
filename={filename.name};
NumberTrajectories = length(filename);

for i=1:1:length(filename)
    data=xlsread(filename{i});
      time=data(:,1)*0.05;
      X=data(:,2:4)/1000;
%     time=data(:,1);
%     X=data(:,2:3);  
    tracks{i} = [time X];
    
end

% Initialize MSDanalyzer
ma = msdanalyzer(3, 'um', 's');
ma = ma.addAll(tracks);
ma = ma.computeMSD;
ma.msd;
MSDdata = ma.getMeanMSD;

%% 
% Plot MSd curves
figure
ma.plotMSD;
legend;
%%
Allmsd = ma.msd;
%% velocity calculation

rows=numel(Allmsd);
velocity=zeros(rows,1);
tau_r=0.04; % 200 nm in 80% gly: 0.37; 200 nm in 80% gly: 0.04
D=0.57;
slopes=zeros(rows,1);
%std=zeros(rows,1);
for i=1:numel(Allmsd)
        a=Allmsd{i};
        t = a(301:500,1);
        x = a(301:500,2);
        [P,S]= polyfit(t,x,1);
        %[fit,delta] = polyval(P,t,S);
        k=P(1);
        slopes(i)=k;
        %std(i)=delta(1);%find how to calculate r-square COD       
        if k-4*D>0
        V=sqrt((k-4*D)./tau_r);
        velocity(i) =V.';       
       end
              
        %plot(t,x,'bo')
        %hold on
        %plot(t,fit,'r-')  
end

