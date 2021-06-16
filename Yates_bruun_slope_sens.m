% Yates model (modified Vitousek) + Bruun Rule (Wolinsky modification)

close all; clc; clear all;
addpath(genpath('/home/isabel/Escritorio/PhD_UoA/Shoreshop_data/Uncertainty_shorelinechange_Vitousek/'));
addpath(genpath('/home/isabel/Escritorio/PhD_UoA'));

%% SLR IPCC 2020-2100 RC P8.5
load('SLR_past_corrected.mat');
time_hist= time;
SLR_hist=SLR;

load('SLR_RCP85_corrected.mat');
%From 2020 to 2100
time_RCP85= time_RCP85(2:10);
SLR_RCP85= SLR_RCP85(2:10);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load Wave hindcast data
load('Wave_hindcast_corrected.mat');
t=hindcast.time;
Hs=hindcast.Hs;

%% shoreline observations to load
load('Shorecast_complete.mat')

tobs=Shore.time;
Yobs=Shore.average-nanmean(Shore.average);

% Cleaning
Yobs(tobs<min(t))=[];
tobs(tobs<min(t))=[];

Yobs(tobs>max(t))=[];
tobs(tobs>max(t))=[];

tobs(isnan(Yobs))=[];
Yobs(isnan(Yobs))=[];

%% Extend wave climate to 2100
t2 = datetime(2100,1,1,0,0,0);
t3 = datestr(t(1)):hours(3):t2;

t3= datenum(t3)';

Hs3=[Hs;Hs;Hs;Hs]; Hs3= Hs3(1:length(t3));
t=t3; Hs=Hs3;

%% Linear fits 2020-2060 %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p_1980=polyfit(time_hist(1:10),SLR_hist(1:10),1);

p_2030= polyfit(time_RCP85(1:2),SLR_RCP85(1:2),1);

p_2060= polyfit(time_RCP85(4:5),SLR_RCP85(4:5),1);

p_2100= polyfit(time_RCP85(8:9),SLR_RCP85(8:9),1);

%% Set sea level rise change vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%RCP8.5 sea level rise%%%%%%%%%%%%%%%%%%%%%%%%%%%
islope= p_1980(1)/1000/8; % slope in meters per 3 hours (days relative)
fslope= p_2100(1) /1000/8 ; % slope in meters per 3 hours (days relative)

dsl= linspace(islope,fslope,length(t)-1)';
%% Set model step and parameters (Yates)

dt=nanmean(diff(t));  % model time step
Hs_bar=nanmean(Hs(:));            % mean wave height

DT=28;   % model time scale
DY=10;   % model shoreline excursion parameter

Nsteps=length(t);
%% Set Bruun parameters 
%tanb = 0.13; %slope at Tairua (Blossier et al . (2017))
tanb = 0.0011; %slope at Miami or similar Anthanasiou (2019) <0.001
tanb2= 0.0013;
tanb3= 0.0016;
%sea level rise rate 3.0 ± 0.4 millimetres  per year for the period 1993–2017
%(Nerem et al. 2018)

%IPCC RCP8.5 by 2100 
%12 mm per year
%S2100= 0.012 / 365/ 8;

%% RUN FORWARD YATES + BRUUN MODEL

Y=NaN(Nsteps,1); Y(1,:)=0;     
Y_b=NaN(Nsteps,1); Y_b(1,:)=0;
Y_b2=NaN(Nsteps,1); Y_b2(1,:)=0;
Y_b3= NaN(Nsteps,1); Y_b3(1,:)=0;
   
for n=1:Nsteps-1
    
    % Y at equilibrium 
    Yeq=-DY*(Hs(n,:).^2-Hs_bar^2)./Hs_bar^2;
    tau=DT*(Hs_bar./Hs(n,:));
    
    % Yates (+ Bruun) wrong formulation
%     Y_b(n+1,:)=Y_b(n,:)+dt./tau.*(Yeq-Y_b(n,:))-(S2100/tanb);
        
    % Yates (unmodified)
    Y(n+1,:)=Y(n,:)+dt./tau.*(Yeq-Y(n,:));
    
    % Yates (+ Bruun), sea level rise change
    Y_b(n+1,:)=Y_b(n,:)+dt./tau.*(Yeq-Y_b(n,:))-(dsl(n)./tanb);
    
    % Yates (+ Bruun), sea level rise change
    Y_b2(n+1,:)=Y_b2(n,:)+dt./tau.*(Yeq-Y_b2(n,:))-(dsl(n)./tanb2);
    
    % Yates (+ Bruun), sea level rise change
    Y_b3(n+1,:)=Y_b3(n,:)+dt./tau.*(Yeq-Y_b3(n,:))-(dsl(n)./tanb3);
end


%% Plots Beach Slope sensibility
figure(3);
sgtitle('Beach slope sensibility')
% subplot(2,1,1);
% plot(tobs(1:10:end),Yobs(1:10:end),'-k'); axis([min(t) max(t) -25 20]); 
% datetick('x','keeplimits');
% ytxt = {'$Y_{obs} [m]$'};
% ylabel(ytxt,'interpreter','latex','FontSize',15);

% subplot(3,1,2);
% %plot monthly 
% plot(t(1:8*15:end),Y(1:8*15:end)); 
% ytxt = {'$Y_{Yates} [m]$'};
% ylabel(ytxt,'interpreter','latex','FontSize',15);
% axis([min(t) max(t) -60 10]); datetick('x','keeplimits');

frec= 8*2000; %100 days

% subplot(2,1,2);
hold on
plot(t,movmean(Y_b,frec),'color',rgb('Violet')); 
plot(t,movmean(Y_b2,frec),'color',rgb('MediumSeaGreen'));
plot(t,movmean(Y_b3,frec),'color',rgb('SteelBlue'));
plot(t,movmean(Y,frec),'color',rgb('Black')); 
ytxt = {'$Y_{Yates+Bruun} [m]$'};
ylabel(ytxt,'interpreter','latex','FontSize',15);
axis([min(t) max(t) -35 -2]); datetick('x','keeplimits');
legend(sprintf('tan = %.4f', tanb),sprintf('tan = %.4f', tanb2),sprintf('tan = %.4f', tanb3),'Yates','Location','southwest')
