# Multi-scale host-aware modeling for analysis and tuning of synthetic gene circuits for bioproduction

PhD dissertation, Fernando N. Santos-Navarro, March 2022.

>“As a physicist, I was used to studying matter that obeys precise mathematical laws. 
But cells are matter that dances.”—Uri Alon, *An Introduction to Systems Biology*

Most of the *OneModel* and *Matlab* code that we have used to generate the simulations and figures in this Thesis is available in this repository.

**Thesis PDF:** **TODO**

**License**: **TODO**

**Contact:** Fernando N. Santos Navarro (<fersann1@upv.es>) and J. Picó (<jpico@ai2.upv.es>)

Please cite using the following BibTex entry: **TODO**

## Abstract

This Thesis was devoted to the multi-scale host-aware analysis and tuning of synthetic gene circuits for bioproduction.
The main objectives were:
  * The development of a reduced-size host-aware model for simulation and analysis purposes.
  * The development of a software toolbox for modeling and simulation, oriented to synthetic biology.
  * The implementation of a multi-scale model that considers the scales relevant to bioproduction (bioreactor, cell, and synthetic circuit).
  * The host-aware analysis of the antithetic controller, as an example of application of the developed tools.
  * The development and experimental validation of robust control laws for continuous bioreactors.

The work presented in this Thesis covers the three scales of the bioproduction process.
The first scale is the bioreactor: this scale considers the macroscopic substrate and biomass dynamics and how these dynamics connect to the internal state of the cells.
The second scale is the host cell: this scale considers the internal dynamics of the cell and the competition for limited shared resources for protein expression.
The third scale is the synthetic genetic circuit: this scale considers the dynamics of expressing exogenous synthetic circuits and the burden they induce on to the host cell.
Finally, as a <<fourth>> scale, part of the Thesis was devoted to developing software tools for modeling and simulation.

This document is divided into seven chapters.
Chapter 1 is an overall introduction to the Thesis work and its justification; it also presents a visual map of the Thesis and lists the main contributions.
Chapter 2 shows the development of the host-aware model (Chapters 4 and 5 make use of this model for their simulations).
Chapter 3 presents \textit{OneModel}: a software tool developed in the Thesis that facilitates modeling and simulation for synthetic biology---in particular, it facilitates the use of the host-aware model---.
Chapter 4 uses the host-aware model to assemble the multi-scale model considering the bioreactor and analyzes the titer, productivity (rate), and yield in expressing an exogenous protein.
Chapter 5 analyzes a more complex circuit, the recently proposed and highly cited antithetic biomolecular controller, using the host-aware model.
Chapter 6 shows the design of nonlinear control strategies that allow controlling the concentration of biomass in a continuous bioreactor in a robust way.
Chapter 7 summarizes and presents the main conclusions of the Thesis.
Appendix 1 shows the theoretical development of the host-aware model.

This Thesis emphasizes the importance of studying cell burden in biological systems since these effects are very noticeable and generate interactions between seemingly unconnected circuits.
The Thesis provides tools to model, simulate and design synthetic genetic circuits taking into account these burden effects and allowing the development of models that connect phenomena in synthetic genetic circuits, ranging from the intracelullar dynamics of gene expression to the macroscopic dynamics of the population of cells inside the bioreactor.
