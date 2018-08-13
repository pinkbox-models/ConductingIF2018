%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% testaxon.m 
% --- for simulating single-compartment membrane response and 
%     multi-compartment axonal spike conduction using the bounded 
%     integrate-and-fire model (compared with the Wang-Buzsaki model)
% Created: Mar 21, 2018 by Go Ashida 
% Revised: Jul 15, 2018 by Go Ashida 
% Revised: Aug 13, 2018 by Go Ashida 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% common time step 
dt = 0.004; % [us]

%% Fig 1: single compartment models (WB/bEIF/sEIF)
% time vector
Tleng = 25; % simulation time length [ms]
TstimS = 10; % stimulus starting time [ms]
TstimE = 20; % stimulus ending time [ms]
t = 0:dt:Tleng; % time vector 

% input vector
iamp = 30; % [pA]
Iinj = zeros(1, length(t));
Iinj((t>=TstimS)) = iamp; % [pA]
Iinj((t>=TstimE)) = 0; 

% calling models
[vW,vB,vS] = singlecompartment(Iinj,dt);

% plotting
figure(1); 
set(gcf,'Position', [100, 500, 560, 420]);
subplot(1,1,1);
cla; hold on; 
plot(t-TstimS,vW,'k-');
plot(t-TstimS,vB,'b-');
plot(t-TstimS,vS,'g-');
xlim([-3,13]);
ylim([-77,27]);
ylabel('potential [mV]');
xlabel('time [ms]');
legend('WB','bEIF','sEIF');
title('single compartment models');
drawnow; 

%% Fig 2: myelinated axon (intracellular current injection)
% morphological parameters
Ndef = 141; % number of nodes
Ddef = 2.0; % [um] axonal diameter 
Ldef = 2.0; % [um] nodal length 
Idef = 200.0; % [um] internodal length 

% time vector
Tleng = 25; % simulation time length [ms]
TstimS = 10; % stimulus starting time [ms]
TstimE = 11; % stimulus ending time [ms]
t = 0:dt:Tleng; % time vector 

% input vector
iamp = 100; % [pA]
Iinj = zeros(Ndef, length(t));
Iinj(20,(t>=TstimS)) = iamp; % [pA]
Iinj(20,(t>=TstimE)) = 0; 

% calling the models
vWB = axonWBintra(Iinj,dt,Ndef,Ddef,Ldef,Idef);
vIF = axonIFintra(Iinj,dt,Ndef,Ddef,Ldef,Idef);

% calculating spike conduction velocity
[m,i1] = max(vWB(40,:));
[m,i2] = max(vWB(90,:));
[m,i3] = max(vIF(40,:));
[m,i4] = max(vIF(90,:));
vel_mWB = ((Ldef+Idef)*50/1000)/(t(i2)-t(i1));
vel_mIF = ((Ldef+Idef)*50/1000)/(t(i4)-t(i3));

% plotting
figure(2); 
set(gcf,'Position', [200, 400, 560*1.5, 420]);

subplot(2,1,1);
cla; hold on; 
plot(t-TstimS,vWB(20,:),'-','color',[0.5,0.5,0.5]); 
plot(t-TstimS,vWB(30,:)-10,'r-'); 
plot(t-TstimS,vWB(40,:)-20,'y-'); 
plot(t-TstimS,vWB(50,:)-30,'g-'); 
plot(t-TstimS,vWB(60,:)-40,'c-'); 
plot(t-TstimS,vWB(70,:)-50,'b-'); 
plot(t-TstimS,vWB(80,:)-60,'m-'); 
plot(t-TstimS,vWB(90,:)-70,'k-'); 
xlim([-1,+6]);
ylim([-160,40]);
ylabel('potential [mV]');
title(sprintf('WB model (myelinated axon): %.2f [m/s]',vel_mWB));

subplot(2,1,2);
cla; hold on; 
plot(t-TstimS,vIF(20,:),'-','color',[0.5,0.5,0.5]); 
plot(t-TstimS,vIF(30,:)-10,'r-'); 
plot(t-TstimS,vIF(40,:)-20,'y-'); 
plot(t-TstimS,vIF(50,:)-30,'g-'); 
plot(t-TstimS,vIF(60,:)-40,'c-'); 
plot(t-TstimS,vIF(70,:)-50,'b-'); 
plot(t-TstimS,vIF(80,:)-60,'m-'); 
plot(t-TstimS,vIF(90,:)-70,'k-'); 
xlim([-1,+6]);
ylim([-160,40]);
ylabel('potential [mV]');
xlabel('time [ms]');
title(sprintf('bEIF model (myelinated axon): %.2f [m/s]',vel_mIF));

drawnow;

%% Fig 3: unmyelinated axon (intracellular current injection)
% morphological parameters
Ndef = 301; % number of nodes
Ddef = 10.0; % [um] axonal diameter 
Ldef = 20.0; % [um] nodal length 
Idef = 0.0; % [um] internodal length (zero for unmyelinated axon)

% time vector
Tleng = 25; % simulation time length [ms]
TstimS = 10; % stimulus starting time [ms]
TstimE = 11; % stimulus ending time [ms]
t = 0:dt:Tleng; % time vector 

% input vector
iamp = 10000; % [pA]
Iinj = zeros(Ndef, length(t));
Iinj(50,(t>=TstimS)) = iamp; % [pA]
Iinj(50,(t>=TstimE)) = 0; 

% calling the models
uWB = axonWBintra(Iinj,dt,Ndef,Ddef,Ldef,Idef);
uIF = axonIFintra(Iinj,dt,Ndef,Ddef,Ldef,Idef);

% calculating spike conduction velocity
[m,i1] = max(uWB(100,:));
[m,i2] = max(uWB(200,:));
[m,i3] = max(uIF(100,:));
[m,i4] = max(uIF(200,:));
vel_uWB = (Ldef*100/1000)/(t(i2)-t(i1));
vel_uIF = (Ldef*100/1000)/(t(i4)-t(i3));

% plotting
figure(3); 
set(gcf,'Position', [300, 300, 560*1.5, 420]);

subplot(2,1,1);
cla; hold on; 
plot(t-TstimS,uWB(50,:),'-','color',[0.5,0.5,0.5]);  
plot(t-TstimS,uWB(75,:)-10,'r-');  
plot(t-TstimS,uWB(100,:)-20,'g-');  
plot(t-TstimS,uWB(125,:)-30,'c-'); 
plot(t-TstimS,uWB(150,:)-40,'b-'); 
plot(t-TstimS,uWB(175,:)-50,'m-'); 
plot(t-TstimS,uWB(200,:)-60,'k-'); 
xlim([-1,+6]);
ylim([-160,40]);
ylabel('potential [mV]');
title(sprintf('WB model (unmyelinated axon): %.2f [m/s]',vel_uWB));

subplot(2,1,2);
cla; hold on; 
plot(t-TstimS,uIF(50,:),'-','color',[0.5,0.5,0.5]);  
plot(t-TstimS,uIF(75,:)-10,'r-');  
plot(t-TstimS,uIF(100,:)-20,'g-');  
plot(t-TstimS,uIF(125,:)-30,'c-'); 
plot(t-TstimS,uIF(150,:)-40,'b-'); 
plot(t-TstimS,uIF(175,:)-50,'m-'); 
plot(t-TstimS,uIF(200,:)-60,'k-'); 
xlim([-1,+6]);
ylim([-160,40]);
ylabel('potential [mV]');
xlabel('time [ms]');
title(sprintf('bEIF model (unmyelinated axon): %.2f [m/s]',vel_uIF));

drawnow;

%% Fig 4: myelinated axon (extracellular current injection)
% morphological parameters
Ndef = 141; % number of nodes
Ddef = 2.0; % [um] axonal diameter 
Ldef = 2.0; % [um] nodal length 
Idef = 200.0; % [um] internodal length 

% time vector
Tleng = 25; % simulation time length [ms]
TstimS = 10; % stimulus starting time [ms]
TstimE = 10.1; % stimulus ending time [ms]
t = 0:dt:Tleng; % time vector 

% input vector
iamp = -1000; % [uA] <- not [pA] 
Iinj = zeros(1, length(t));
Iinj(t>=TstimS) = iamp; % [pA]
Iinj(t>=TstimE) = 0; 

% distances between current source and nodes
% assuming that the current source is located at 1 mm away from node #20
rz = 1.0; % [mm]
rx = ( (1:Ndef)-20 )' * (Ldef+Idef) /1000; % [mm]
rc = sqrt( rz.*rz + rx.*rx );

% calling the models
eWB = axonWBextra(Iinj,rc,dt,Ndef,Ddef,Ldef,Idef);
eIF = axonIFextra(Iinj,rc,dt,Ndef,Ddef,Ldef,Idef);

% calculating spike conduction velocity
[m,i1] = max(eWB(40,:));
[m,i2] = max(eWB(90,:));
[m,i3] = max(eIF(40,:));
[m,i4] = max(eIF(90,:));
vel_eWB = ((Ldef+Idef)*50/1000)/(t(i2)-t(i1));
vel_eIF = ((Ldef+Idef)*50/1000)/(t(i4)-t(i3));

% plotting
figure(4); 
set(gcf,'Position', [400, 200, 560*1.5, 420]);

subplot(2,1,1);
cla; hold on; 
plot(t-TstimS,eWB(20,:),'-','color',[0.5,0.5,0.5]); 
plot(t-TstimS,eWB(30,:)-10,'r-'); 
plot(t-TstimS,eWB(40,:)-20,'y-'); 
plot(t-TstimS,eWB(50,:)-30,'g-'); 
plot(t-TstimS,eWB(60,:)-40,'c-'); 
plot(t-TstimS,eWB(70,:)-50,'b-'); 
plot(t-TstimS,eWB(80,:)-60,'m-'); 
plot(t-TstimS,eWB(90,:)-70,'k-'); 
xlim([-1,+6]);
ylim([-160,40]);
ylabel('potential [mV]');
title(sprintf('WB model (extracellular stimulation): %.2f [m/s]',vel_eWB));

subplot(2,1,2);
cla; hold on; 
plot(t-TstimS,eIF(20,:),'-','color',[0.5,0.5,0.5]); 
plot(t-TstimS,eIF(30,:)-10,'r-'); 
plot(t-TstimS,eIF(40,:)-20,'y-'); 
plot(t-TstimS,eIF(50,:)-30,'g-'); 
plot(t-TstimS,eIF(60,:)-40,'c-'); 
plot(t-TstimS,eIF(70,:)-50,'b-'); 
plot(t-TstimS,eIF(80,:)-60,'m-'); 
plot(t-TstimS,eIF(90,:)-70,'k-'); 
xlim([-1,+6]);
ylim([-160,40]);
ylabel('potential [mV]');
xlabel('time [ms]');
title(sprintf('bEIF model (extracellular stimulation): %.2f [m/s]',vel_eIF));

drawnow;

%% Fig 5 (part 1): single compartment auditory nerve
% time vector
Tleng = 25; % simulation time length [ms]
TstimS = 10; % stimulus starting time [ms]
TstimE = 20.5; % stimulus ending time [ms]
t = 0:dt:Tleng; % time vector 

% input vector
iamp = 40; % [pA]
iext = zeros(1, length(t));
iext((t>=TstimS)) = iamp; % [pA]
iext((t>=TstimE)) = 0; 

% calling models
vL = ANsingle(iext,dt,0); % low freq AN
vH = ANsingle(iext*2,dt,1); % high freq AN

% plotting
figure(5); 
set(gcf,'Position', [600, 300, 560*1.5, 420*1.5]);
subplot(3,2,1);
cla; hold on; 
plot(t-TstimS,vL,'r-');
plot(t-TstimS,vH,'b-');
xlim([-3,16]);
ylim([-77,35]);
ylabel('potential [mV]');
xlabel('time [ms]');
title('single compartment AN');
legend('low freq','high freq','Location','NorthWest');

subplot(3,2,2);
cla; hold on; 
plot((t-TstimS)*0.75,vL,'r-');
plot((t-TstimS)*1.50,vH,'b-');
set(gca,'xtick',0);
xlim([-3,16]);
ylim([-77,35]);
title('single compartment AN (time scaled)');
legend('low freq','high freq','Location','NorthWest');

drawnow;

%% Fig 5 (part 2): multi-compartment auditory nerve axon  
% morphological parameter (other parameters are defined in axonAN.m)
Ndef = 40; % number of nodes

% time vector
Tleng = 25; % simulation time length [ms]
TstimS = 10; % stimulus starting time [ms]
TstimE = 11; % stimulus ending time [ms]
t = 0:dt:Tleng; % time vector 

% input vector
iamp = 60; % [pA]
iext = zeros(Ndef, length(t));
iext(1,(t>=TstimS)) = iamp; % [pA]
iext(1,(t>=TstimE)) = 0; 

% calling the models
vL = ANaxon(iext,dt,Ndef,0); % low freq model 
vH = ANaxon(iext,dt,Ndef,1); % high freq model 

% spike conduction velocity
[m,i1] = max(vL(10,:));
[m,i2] = max(vL(30,:));
[m,i3] = max(vH(10,:));
[m,i4] = max(vH(30,:));
Ldef = 2.0; % nodal length 
Idef0 = 350; % internodal length (low freq)
Idef1 = 450; % internodal length (high freq)
velL = ((Ldef+Idef0)*20/1000)/(t(i2)-t(i1));
velH = ((Ldef+Idef1)*20/1000)/(t(i4)-t(i3));

% plotting
figure(5); 
set(gcf,'Position', [600, 300, 560*1.5, 420*1.5]);

subplot(3,2,[3,4]);
cla; hold on; 
plot(t-TstimS,vL(1,:),'-','color',[0.5,0.5,0.5]); 
plot(t-TstimS,vL(6,:)-10,'r-'); 
plot(t-TstimS,vL(11,:)-20,'y-'); 
plot(t-TstimS,vL(16,:)-30,'g-'); 
plot(t-TstimS,vL(21,:)-40,'c-'); 
plot(t-TstimS,vL(26,:)-50,'b-'); 
plot(t-TstimS,vL(31,:)-60,'m-'); 
plot(t-TstimS,vL(36,:)-70,'k-'); 
xlim([-1,+4]);
ylim([-160,40]);
ylabel('potential [mV]');
title(sprintf('low frequency AN axon model: %.2f [m/s]',velL));

subplot(3,2,[5,6]);
cla; hold on; 
plot(t-TstimS,vH(1,:),'-','color',[0.5,0.5,0.5],'LineWidth',1); 
plot(t-TstimS,vH(6,:)-10,'r-'); 
plot(t-TstimS,vH(11,:)-20,'y-'); 
plot(t-TstimS,vH(16,:)-30,'g-'); 
plot(t-TstimS,vH(21,:)-40,'c-'); 
plot(t-TstimS,vH(26,:)-50,'b-'); 
plot(t-TstimS,vH(31,:)-60,'m-'); 
plot(t-TstimS,vH(36,:)-70,'k-'); 
xlim([-1,+4]);
ylim([-160,40]);
ylabel('potential [mV]');
xlabel('time [ms]');
title(sprintf('high frequency AN axon model: %.2f [m/s]',velH));

drawnow;

