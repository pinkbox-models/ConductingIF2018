function [vW,vB,vS] = singlecompartment(Iinj,dt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Single compartment models 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input 
%  Iinj : injected current [pA] (vector)
%  dt   : time step [ms] 
% Output 
%  vW : membrane potential (Wang-Buzsaki model)
%  vB : membrane potential (bounded exponential integrate-and-fire model)
%  vS : membrane potential (standard exponential integrate-and-fire model)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes
% + See the references below for the details of the models.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References 
% Ashida G, Nogueira W (2018) 
%  "Spike-conducting integrate-and-fire model" 
%  eNeuro (to be published online)
% Fourcaud-Trocme N et al. (2003)
%  "How spike generation mechanisms determine the neuronal response to 
%   fluctuating inputs", J Neurosci 23: 11628-11640 
% Wang XJ, Buzsaki G (1996) 
%  "Gamma oscillation by synaptic inhibition in a hippocampal 
%   interneuronal network model", J Neurosci 16: 6402-6413 
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

%% common membrane parameters
Smemb = 1000; % [um^2] surface area of the membrane 
Cmemb = 1.0;  % [uF/cm^2] membrane capacitance density 
Cm = Cmemb * Smemb * 1e-8; % [uF] membrane capacitance 

%% model-specific parameters
% WB model
pW = struct();
pW.GL = 0.1; % [mS/cm^2] leak conductance density 
pW.GN = 35.0; % [mS/cm^2] Na conductance density
pW.GK = 15.0; % [mS/cm^2] K conductance density 
pW.EN = +55.0; % [mV] Na reversal potential 
pW.EK = -90.0; % [mV] K reversal potential 
pW.EL = -65.0; % [mV] leak reversal potential 
pW.gN = pW.GN * Smemb * 1e-8; % Na conductance [mS]
pW.gK = pW.GK * Smemb * 1e-8; % K conductance [mS]
pW.gL = pW.GL * Smemb * 1e-8; % leak conductance [mS]

% bounded EIF model
pB = struct();
pB.GL = 0.1; % [mS/cm^2] leak conductance density 
pB.EL = -65.3; % [mV] leak reversal potential 
pB.Vth = -60.2; % [mV] spike threshold 
pB.Kth = 3.5; % [mV] slope factor 
pB.Ath = 520; % [no unit] (Note: max(g_depolarize) = gL*Ath)
pB.Vsp = +10; % [mV] spike-detecting threshold 
pB.spT = 0.60; % [ms] repolarization current time constant
pB.spA = 90; % [no unit] (Note: max(g_repolarize)= gL*spA)
pB.Vre = pB.EL + 20; % [mV] voltage for assuring repolarization
pB.gL = pB.GL * Smemb * 1e-8; % leak conductance [mS]

% standard EIF model
pS = struct();
pS.GL = 0.1; % [mS/cm^2] leak conductance density 
pS.EL = -65.3; % [mV] leak reversal potential 
pS.Vth = -60.2; % [mV] spike threshold 
pS.Kth = 3.5; % [mV] slope factor 
pS.Vsp = +15; % [mV] spike-detecting threshold 
pS.Tref = 2.8; % [ms] refractory period 
pS.gL = pS.GL * Smemb * 1e-8; % leak conductance [mS]

%% vectors for storing variables
% WB model
vW = zeros(1, length(Iinj)); % [mV] membrane potential 
mW = zeros(1, length(Iinj)); % Na activation variable
hW = zeros(1, length(Iinj)); % Na inactivation variable
nW = zeros(1, length(Iinj)); % K activation variable 
inW = zeros(1, length(Iinj)); % Na current
ikW = zeros(1, length(Iinj)); % K current
ilW = zeros(1, length(Iinj)); % leak current 

% bEIF model
vB = zeros(1,length(Iinj)); % [mV] membrane potential 
ilB = zeros(1, length(Iinj)); % leak current 
idB = zeros(1, length(Iinj)); % depolarizing current
irB = zeros(1, length(Iinj)); % repolarizing current 
xrB = zeros(1, length(Iinj)); % repolarizing conductance (main)
yrB = zeros(1, length(Iinj)); % repolarizing conductance (sub)
spReadyB = 1; % ready-to-spike flag
spStartB = 0; % start-of-spike flag

% sEIF model
vS = zeros(1,length(Iinj)); % [mV] membrane potential 
ilS = zeros(1, length(Iinj)); % leak current 
idS = zeros(1, length(Iinj)); % leak current 
NrefS = round(pS.Tref/dt); % refractory period in steps 
rCountS = 0; % counter for refractory period 

% parameters for calculating alpha functions 
pB.ee = exp(-dt/pB.spT); % decay factor 
pB.aa = pB.spA * exp(1.0) / pB.spT; % amplitude factor

% initial values
vW(:,1) = pW.EL;  % initial membrane potential (WB model) 
mW(:,1) = WBalphaM(pW.EL) / (WBalphaM(pW.EL) + WBbetaM(pW.EL)) ; % initial m
hW(:,1) = WBalphaH(pW.EL) / (WBalphaH(pW.EL) + WBbetaH(pW.EL)) ; % initial h
nW(:,1) = WBalphaN(pW.EL) / (WBalphaN(pW.EL) + WBbetaN(pW.EL)) ; % initial n
vB(:,1) = pB.EL;  % initial membrane potential (bEIF model) 
vS(:,1) = pS.EL;  % initial membrane potential (sEIF model) 

%% calculate membrane response step-by-step 
for j=1:length(Iinj)-1

    % ionic currents (WB model) g[mS] * V[mV] = I[uA]
    inW(j) = pW.gN .* mW(:,j).^3 .* hW(:,j) .* ( pW.EN - vW(j) ); 
    ikW(j) = pW.gK .* nW(:,j).^4            .* ( pW.EK - vW(j) ); 
    ilW(j) = pW.gL                          .* ( pW.EL - vW(j) ); 

    % ionic currents (bEIF model) 
    ilB(j) = pB.gL .* ( pB.EL - vB(j) ); 
    idB(j) = pB.gL .* pB.Kth  .* pB.Ath ./ (1 + pB.Ath * exp( -(vB(j)-pB.Vth)/pB.Kth )); 
    irB(j) = pB.gL .* xrB(j) .* ( pB.EL - vB(j) );

    % ionic currents (sEIF model)
    ilS(j) = pS.gL * ( pS.EL-vS(j) ); 
    idS(j) = pS.gL * pS.Kth * exp( (vS(j)-pS.Vth)/pS.Kth ); 

    % derivatives I[uA] / C[uF] * dt[ms] = dv[mV]    
    dvW_dt = ( inW(j) + ikW(j) + ilW(j) + Iinj(j)*1e-6 ) ./ Cm;  
    dmW_dt = (1-mW(j)).* WBalphaM(vW(j)) - mW(j).*WBbetaM(vW(j));
    dhW_dt = (1-hW(j)).* WBalphaH(vW(j)) - hW(j).*WBbetaH(vW(j));
    dnW_dt = (1-nW(j)).* WBalphaN(vW(j)) - nW(j).*WBbetaN(vW(j));
    dvB_dt = ( ilB(j) + idB(j) + irB(j) + Iinj(j)*1e-6 ) ./ Cm;  
    dvS_dt = ( ilS(j) + idS(j) + Iinj(j)*1e-6 ) ./ Cm;  

    % calculate next step 
    vW(j+1) = dvW_dt * dt + vW(j); 
    mW(:,j+1) = mW(:,j) + dmW_dt * dt; 
    hW(:,j+1) = hW(:,j) + dhW_dt * dt; 
    nW(:,j+1) = nW(:,j) + dnW_dt * dt; 
    vB(j+1) = dvB_dt * dt + vB(j); 
    vS(j+1) = dvS_dt * dt + vS(j); 

    % check for threshold crossing (bEIF model)
    if(spReadyB==1) % if ready to spike, then check for spike occurrence
      if(vB(j+1)>=pB.Vsp) 
        spStartB = 1; % now started spiking
        spReadyB = 0; % thus not ready for spiking
      else 
        spStartB = 0; % not started spiking yet
        spReadyB = 1; % thus still ready for spiking
      end
    else % if in repolarizing phase, then check whether voltage is back near rest
      if(vB(j+1)<=pB.Vre)
        spStartB = 0; % not starting spiking, of course 
        spReadyB = 1; % and now ready for next spike
      else
        spStartB = 0; % not starting spiking, of course 
        spReadyB = 0; % not yet ready for next spike
      end
    end

    % step for alpha function (bEIF model)
    yrB(j+1) = pB.ee * yrB(j) + pB.aa .* spStartB; 
    xrB(j+1) = pB.ee * xrB(j) + dt * yrB(j+1);

    % check for threshold crossing (of sEIF model)
    if(rCountS>0); 
      rCountS = rCountS-1;  
      vS(j+1) = pS.EL;
    elseif(vS(j+1)>=pS.Vsp) 
      rCountS = NrefS; 
      vS(j+1) = pS.Vsp;
    end 

end 

%% WB activation/inactivation functions
function x = WBalphaM(v)
  x = 5 * 0.1 * (v+35) ./ ( 1 - exp(-(v+35)/10) ); 
function x = WBalphaH(v)
  x = 0.35 * exp(-(v+58)/20);
function x = WBalphaN(v)
  x = 0.05 * (v+34) ./ ( 1 - exp(-(v+34)/10) );
function x = WBbetaM(v)
  x = 5 * 4.0 * exp(-(v+60)/18);
function x = WBbetaH(v)
  x = 5.0 ./ ( 1 + exp(-(v+28)/10) ); 
function x = WBbetaN(v)
  x = 0.625 * exp(-(v+44)/80);
