# RichnessPatternsAlongEnvironmentalGradients
The model presented in the main text of the Journal of Biogeography paper “The origin of species richness patterns along environmental gradients: uniting explanations based on time, diversification rate and carrying capacity" by Mikael Pontarp and John J. Wiens. The model is implemented in its basic form as MATLAB code. To run the model “main_simulation_prog_5pach.m” should be executed in the same directory as the other m-files. The code is commented within each m-file, below we present general features and key components of the implementation.

# Model code, general instructions
This is the main function that is executed to run the model. Model parameters and initial conditions are initiated at lines 36-81. The function requires five input arguments, each reflecting a given parameter (see lines 3-7). Such initiation also allows for the possibility of atomized model output for different parts of parameter space or different simulation realizations (e.g. replicates). 

For each set of simulations (defined as a combination of biotic and abiotic conditions, specified below), the code simulates alternating phases of individual reproduction and dispersal for 100,000 generations (time-step) and each simulation can be replicated (modifications to code required). At the beginning of each simulation realization, a habitat at the extreme end of the gradient is seeded with 10 monomorphic individuals. During reproduction, each individual reproduced according to its fitness and each offspring inherited the same trait values as their parent (asexual reproduction) unless the offspring mutated (source of phenotypic variation). All offspring are born into the habitat of their parent but dispersed with a probability (d) during the dispersal phase to one of the neighboring habitats according to a stepping-stone dispersal algorithm. The simulation follows all individuals, calculates their fitness, and allows them to reproduce and disperse. As a result, the simulation output is a distribution of individuals in trait space and geographical space for each time step. The mutation process (with offspring values close to parental values), together with the evolutionary process (driven by the fitness-generating function), generates a clustered distribution of trait values along the trait axes (Fig. 2). We treat these
clusters of similar individuals as species (Pontarp et al.,
2012, 2015). See also Appendix S1 of the published paper for details.

# main_simulation_prog_5patch.m
Depending on input arguments this function runs the simulation as described above.

# fitfunc.m
Called by the main script, used for computing individual fitness. 

# speccheck.m
Called by the main script. Identifies clusters of phenotypically similar individuals and assigns "speciation" events and "species" identity accordingly. A trait-based speciation concept is used.   

# treeconstruc.m
Analyses speciation history and phylogenetic relatedness among "species". A phylogenetic tree is produced. 
