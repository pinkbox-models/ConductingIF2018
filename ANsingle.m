function v = ANsingle(Iinj,dt,F)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Auditory nerve axonal spike conduction model 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input 
%  Iinj : injected current [pA] (vector)
%  dt   : time step [ms] 
%  F    : 0=low frequency model, 1=high frequency model
% Output
%  v    : membrane potential at each node at each time step 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reference 
% Ashida G, Nogueira W (2018) 
%  "Spike-conducting integrate-and-fire model" 
%  eNeuro (to be published online)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Revisions 
% Created (ver. 0.9.0): Mar 21, 2018 by GA
% Revised (ver. 0.9.1): Jul 15, 2018 by GA
% Revised (ver. 0.9.2): Aug 13, 2018 by GA 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you find a bug, please report to GA at go.ashida@uni-oldenburg.de
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% Copyright 2018 Go Ashida (go.ashida@uni-oldenburg.de) %%%%%%%%%%%%%
% Permission is hereby granted under the Apache License, Version 2.0; 
% Users of this file must be in compliance with this license, a copy of 
% which may be obtained at http://www.apache.org/licenses/LICENSE-2.0
% This file is provided on an "AS IS" basis, WITHOUT WARRANTIES OR 
% CONDITIONS OF ANY KIND, either express or implied.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% physiological parameters
% frequency-independent parameters
Smemb = 1000.0; % [um^2]
Cmemb = 1.0; % [uF/cm^2] membrane capacitance density 
EL = -65.3; % [mV] leak reversal potential 
Vth = -50.0; % [mV] spike threshold 
Kth = 3.5; % [mV] slope factor 
Ath = 520; % [no unit] (Note: max(g_depolarize) = gL*Ath)
Vsp = +10; % [mV] spike-detecting threshold 
spT = 0.60; % [ms] repolarization current time constant
spA = 90; % [no unit] (Note: max(g_repolarize)= gL*spA)
Cm = Cmemb * Smemb * 1e-8; % [uF] membrane capacitance 

% frequency-dependent parameter
if F==0 % low freq
  GL = 0.2; % [mS/cm^2] leak conductance density 
elseif F==1 % high freq
  GL = 0.4; % [mS/cm^2] leak conductance density 
end
gL = GL * Smemb * 1e-8; % leak conductance [mS]


%% internal parameters for spike generation and repolarization
Vre = EL + 20; % [mV] voltage for assuring repolarization
spReady = 1; % ready-to-spike flag
spStart = 0; % start-of-spike flag

% parameters for calculating alpha functions 
ee = exp(-dt/spT); % decay factor 
aa = spA * exp(1.0) / spT; % amplitude factor

%% vectors for storing variables
v = zeros(1, length(Iinj)); % [mV] membrane potential 
iL = zeros(1, length(Iinj)); % leak current 
iD = zeros(1, length(Iinj)); % depolarizing current
iR = zeros(1, length(Iinj)); % repolarizing current 
xR = zeros(1, length(Iinj)); % repolarizing conductance (main)
yR = zeros(1, length(Iinj)); % repolarizing conductance (sub)

% initial values
v(1) = EL;  % initial membrane potential

%% calculate membrane response step-by-step 
for j=1:length(Iinj)-1

    % ionic currents g[mS] * V[mV] = I[uA]
    iL(j) = gL .* ( EL - v(j) ); % leak
    iD(j) = gL .* Kth  .* Ath ./ (1 + Ath * exp( -(v(j)-Vth)/Kth )); % depolarization
    iR(j) = gL .* xR(j) .* ( EL - v(j) ); % repolarization

    % derivatives I[uA] / C[uF] * dt[ms] = dv[mV]    
    dv_dt = ( iL(j) + iD(j) + iR(j) + Iinj(j)*1e-6 ) ./ Cm;  

    % calculate next step 
    v(j+1) = dv_dt * dt + v(j); 
    
    % check for threshold crossing (bEIF model)
    if(spReady==1) % if ready to spike, then check for spike occurrence
      if(v(j+1)>=Vsp) 
        spStart = 1; % now started spiking
        spReady = 0; % thus not ready for spiking
      else 
        spStart = 0; % not started spiking yet
        spReady = 1; % thus still ready for spiking
      end
    else % if in repolarizing phase, then check whether voltage is back near rest
      if(v(j+1)<=Vre)
        spStart = 0; % not starting spiking, of course 
        spReady = 1; % and now ready for next spike
      else
        spStart = 0; % not starting spiking, of course 
        spReady = 0; % not yet ready for next spike
      end
    end

    % step for alpha function 
    yR(j+1) = ee * yR(j) + aa .* spStart; 
    xR(j+1) = ee * xR(j) + dt * yR(j+1);

end 
