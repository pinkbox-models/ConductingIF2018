function v = axonWBintra(Iinj,dt,Nc,Dc,Lc,Ic)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Axonal spike conduction simulated with Wang-Buzsaki model 
% stimulated with intracellular current injection
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
% + See the references below for the details of the model.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References 
% Ashida G, Nogueira W (2018) 
%  "Spike-conducting integrate-and-fire model" 
%  eNeuro (to be published online)
% Wang XJ, Buzsaki G (1996) 
%  "Gamma oscillation by synaptic inhibition in a hippocampal 
%   interneuronal network model", J Neurosci 16: 6402-6413 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Revisions 
% Created (ver. 0.9.0): Jun 21, 2018 by GA
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
GN = 35.0; % [mS/cm^2] Na conductance density
GK = 15.0; % [mS/cm^2] K conductance density 
EN = +55.0; % [mV] Na reversal potential 
EK = -90.0; % [mV] K reversal potential 
EL = -65.0; % [mV] leak reversal potential 
Cm = Cmemb * Scomp * 1e-8; % [uF] membrane capacitance 
gN = GN * Scomp * 1e-8; % Na conductance [mS]
gK = GK * Scomp * 1e-8; % K conductance [mS]
gL = GL * Scomp * 1e-8; % leak conductance [mS]
Rax = 1.0; % [MOhm.um] axial resistivity 

%% vectors for storing variables
v = zeros(Nc, length(Iinj)); % [mV] membrane potential 
m = zeros(Nc, length(Iinj)); % Na activation variable
h = zeros(Nc, length(Iinj)); % Na inactivation variable
n = zeros(Nc, length(Iinj)); % K activation variable 
iN = zeros(Nc, length(Iinj)); % Na current
iK = zeros(Nc, length(Iinj)); % K current
iL = zeros(Nc, length(Iinj)); % leak current 

% initial values
v(:,1) = EL;  % initial membrane potential
m(:,1) = WBalphaM(EL) / (WBalphaM(EL) + WBbetaM(EL)) ; % initial m
h(:,1) = WBalphaH(EL) / (WBalphaH(EL) + WBbetaH(EL)) ; % initial h
n(:,1) = WBalphaN(EL) / (WBalphaN(EL) + WBbetaN(EL)) ; % initial n

%% diffusion matrix
[Adiff, Bdiff] = diffusionmatrix(dt, Rax, Nc, Dcomp, Lcomp, Icomp, Cm);

%% calculate membrane response step-by-step 
for j=1:length(Iinj)-1

    % ionic currents  g[mS] * V[mV] = I[uA]
    iN(:,j) = gN .* m(:,j).^3 .* h(:,j) .* ( EN - v(:,j) ); 
    iK(:,j) = gK .* n(:,j).^4           .* ( EK - v(:,j) ); 
    iL(:,j) = gL                        .* ( EL - v(:,j) ); 

    % derivatives I[uA] / C[uF] * dt[ms] = dv[mV]    
    dv_dt = ( iN(:,j) + iK(:,j) + iL(:,j) + Iinj(:,j)*1e-6 ) ./ Cm;  
    dm_dt = (1-m(:,j)).* WBalphaM(v(:,j)) - m(:,j).*WBbetaM(v(:,j));
    dh_dt = (1-h(:,j)).* WBalphaH(v(:,j)) - h(:,j).*WBbetaH(v(:,j));
    dn_dt = (1-n(:,j)).* WBalphaN(v(:,j)) - n(:,j).*WBbetaN(v(:,j));

    % calculate next step 
    v(:,j+1) = Adiff \ ( dv_dt * dt + Bdiff*v(:,j)  ); 
    
    % calculate next step (for WB channels)
    m(:,j+1) = m(:,j) + dm_dt * dt; 
    h(:,j+1) = h(:,j) + dh_dt * dt; 
    n(:,j+1) = n(:,j) + dn_dt * dt; 

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
