------------------------------------------------------------------------------------
  Spike-Conducting Integrate-and-Fire Model -- Matlab Implementation (ver. 0.9.2)
------------------------------------------------------------------------------------

%%% Versions %%% 

- Ver. 0.9.0 (Mar. 21, 2018): Created to be submitted as "Extended Data" to the journal with the manuscript 
- Ver. 0.9.1 (Jul. 15, 2018): Revised with additional code for extracellular stimulation 
- Ver. 0.9.2 (Aug. 13, 2018): Revised with additional comments 


%%% Author %%% 

Go Ashida (University of Oldenburg) go.ashida@uni-oldenburg.de


%%% Abbreviations %%% 

WB : Wang-Buzsaki (model) 
EIF: exponential integrate-and-fire (model) 
bEIF: bounded EIF (model) 
sEIF: standard EIF (model) 
AN: auditory nerve 


%%% Contents %%% 

- testaxon.m : Sample code for plotting single- and multi-compartment spiking responses of the models. 
               This script internally uses the following .m files.

- singlecompartment.m : Single compartment models (WB, bEIF, and sEIF) 

- axonWBintra.m : Intracellularly stimulated multi-compartment spike-conducting axon based on WB model. 
- axonIFintra.m : Intracellularly stimulated multi-compartment spike-conducting axon based on bEIF model.

- axonWBextra.m : Extracellularly stimulated multi-compartment spike-conducting axon based on WB model. 
- axonIFextra.m : Extracellularly stimulated multi-compartment spike-conducting axon based on bEIF model. 

- ANsingle.m : Single compartment bEIF model with AN parameters. 
- ANaxon.m : Multi-compartment spike-conducting bEIF model axon with AN parameters. 

- diffusionmatrix.m : Code for calculating diffusion matrices (internally called by axonXXX.m files) 

+ Notes: See each program file and the reference below for more detailed descriptions. 


%%% Reference %%% 

Ashida G, Nogueira W (2018) 
"Spike-conducting integrate-and-fire model" 
eNeuro (to be published online)

