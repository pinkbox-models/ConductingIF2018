function v = axonIFintra(Iinj,dt,Nc,Dc,Lc,Ic)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Axonal spike conduction simulated with bounded integrate-and-fire (bEIF) 
% model stimulated with intracellular current injection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input 
%  Iinj : injected current to each node at each time step [pA] (2D-vector)
%  dt   : time step [ms] 
%  Nc   : number of compartments 
%  Dc   : axonal diameter [um] 
%  Lc   : nodal length [um] 
%  Ic   : internodal length [um] 
% Output
%  v    : membrane potential at each node at each time step [mV] (2D-vector)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes
% + See the reference below for the details of the model.
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

%% compartment parameters
Dcomp = repmat(Dc, Nc, 1); % [um] axonal diameter
Lcomp = repmat(Lc, Nc, 1); % [um] nodal length 
Icomp = repmat(Ic, Nc-1, 1); % [um] internodal length
Scomp = pi * Dcomp .* Lcomp; % [um2] surface area 

%% physiological parameters
Cmemb = 1.0; % [uF/cm^2] membrane capacitance density 
GL = 0.1; % [mS/cm^2] leak conductance density 
EL = -65.3; % [mV] leak reversal potential 
Vth = -60.2; % [mV] spike threshold 
Kth = 3.5; % [mV] slope factor 
Ath = 520; % [no unit] (Note: max(g_depolarize) = gL*Ath)
Vsp = +10; % [mV] spike-detecting threshold 
spT = 0.60; % [ms] repolarization current time constant
spA = 90; % [no unit] (Note: max(g_repolarize)= gL*spA)
Cm = Cmemb * Scomp * 1e-8; % [uF] membrane capacitance 
gL = GL * Scomp * 1e-8; % leak conductance [mS]
Rax = 1.0; % [MOhm.um] axial resistivity 

%% internal parameters for spike generation and repolarization
Vre = EL + 20; % [mV] voltage for assuring repolarization
spReady = ones(Nc,1); % ready-to-spike flag
spStart = zeros(Nc,1); % start-of-spike flag

% parameters for calculating alpha functions 
ee = exp(-dt/spT); % decay factor 
aa = spA * exp(1.0) / spT; % amplitude factor

%% vectors for storing variables
v = zeros(Nc, length(Iinj)); % [mV] membrane potential 
iL = zeros(Nc, length(Iinj)); % leak current 
iD = zeros(Nc, length(Iinj)); % depolarizing current
iR = zeros(Nc, length(Iinj)); % repolarizing current 
xR = zeros(Nc, length(Iinj)); % repolarizing conductance (main)
yR = zeros(Nc, length(Iinj)); % repolarizing conductance (sub)

% initial values
v(:,1) = EL;  % initial membrane potential

%% diffusion matrix
[Adiff, Bdiff] = diffusionmatrix(dt, Rax, Nc, Dcomp, Lcomp, Icomp, Cm);

%% calculate membrane response step-by-step 
for j=1:length(Iinj)-1

    % ionic currents g[mS] * V[mV] = I[uA]
    iL(:,j) = gL .* ( EL - v(:,j) ); % leak 
    iD(:,j) = gL .* Kth  .* Ath ./ (1 + Ath * exp( -(v(:,j)-Vth)/Kth )); % depolarization
    iR(:,j) = gL .* xR(:,j) .* ( EL - v(:,j) ); % repolarization

    % derivatives I[uA] / C[uF] * dt[ms] = dv[mV]    
    dv_dt = ( iL(:,j) + iD(:,j) + iR(:,j) + Iinj(:,j)*1e-6 ) ./ Cm;  

    % calculate next step 
    v(:,j+1) = Adiff \ ( dv_dt * dt + Bdiff*v(:,j)  ); 
    
    % check for threshold crossing 
    for n=1:Nc
      if(spReady(n)==1) % if ready to spike, then check for spike occurrence
        if(v(n,j+1)>=Vsp) 
          spStart(n) = 1; % now started spiking
          spReady(n) = 0; % thus not ready for spiking
        else 
          spStart(n) = 0; % not started spiking yet
          spReady(n) = 1; % thus still ready for spiking
        end
      else % if in repolarizing phase, then check whether voltage is back near rest
        if(v(n,j+1)<=Vre)
          spStart(n) = 0; % not starting spiking, of course 
          spReady(n) = 1; % and now ready for next spike
        else
          spStart(n) = 0; % not starting spiking, of course 
          spReady(n) = 0; % not yet ready for next spike
        end
      end
    end % end of for loop
 
    % step for alpha function 
    yR(:,j+1) = ee * yR(:,j) + aa .* spStart; 
    xR(:,j+1) = ee * xR(:,j) + dt * yR(:,j+1);

end 
