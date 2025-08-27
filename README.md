The code simulates the neuronal activity of an Leak-integrated and fire excitatory ring network coupled to an astrocyte network following dynamics presented in Menesse et. al (2025). 

# RUNNING MODEL

Example notebook code for building the network with the specific topology used in Menesse et. al (2025) and running a simulation.
The same notebook contained a cell for plot the time series of many variables including a raster plot neuronal activity.

# JULIA FUNCTIONS
The script AstroCircuit_Functions.jl contains functions for setting up the stimuli protocol and solve system of ODE's using a simple euler method.

The output of the script gives:
  - File with spike time data: 
    Columns:
    1- Neuron ID.
    2- Spike time (t_spk).
    3- Fraction of activated neurotransmiter at spike time "y(t_spk)".
    4- Release probabiilty of pre-synaptic component at spike time "U_SE(T_spk)".
    5- Fraction of pre-synaptic neurotransmiter resource available at spike time "x^{pre}(T_spk)".
    6- Membrane potentail of the post-synaptic neuron at spike time "V(T_spk)"
    7- Difference between current time (t_spk) and previous spike.

# Reference:
Menesse, G., Milan, A. P., & Torres, J. J. (2025). Astrocyte-Mediated Higher-Order Control of Synaptic Plasticity. DOI: https://doi.org/10.48550/arXiv.2507.07693
