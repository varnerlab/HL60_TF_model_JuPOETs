## HL60 Gene Regulatory Network Model


### Background ###
The HL60 cell line is an attractive, comprehensive model for understanding how RA-induced differentiation and proliferation mechanisms operate. HL60 cells are bipotent myelomonocttic precursors, capable of being induced to differnetiate into monocytes or granulocytes. Treatments with all-trans retinoic acid (RA) induces differentiaition of HL60 along the granuloctyic linage into neutrophil-like cells. In this study, we developed an ensemble of experimentally validated model for gene expression. We formulated a gene regulatory network that control cell-arrest, and linage commitement. Using the networkd, we combined oridinary differential euations with logical rules to produce a compact yet predictivemodel. To confirm our network structure, we utilized experimental measurements of phenotypic markers that are regulated by one or a combination of the molecules in our gene regulatory network.


### Installation

You can download this repository as a zip file, or clone or pull it by using the command:

	git pull https://github.com/varnerlab/HL60_TF_model_JuPOETs.git

or

	git clone https://github.com/varnerlab/HL60_TF_model_JuPOETs.git

### Model code and parameter ensemble
The HL60 GRN equations were implemented in [Julia](http://julialang.org) and solved using the od23s routine of the [ODE package](https://github.com/JuliaDiffEq/ODE.jl). The model code and parameter ensemble is freely available under an [MIT software license](https://opensource.org/licenses/MIT).

The model equations are encoded in ``Balances.jl`` which is called by the ``SolveBalances.jl`` driver function. The user should not directly call ``SolveBalances.jl``. Rather, multiple parameter sets can be simulated by calling the driver function from a script. The kinetic and other model parameters are encoded in ``DataDictionary.jl`` as a dictionary. The parameters stored in this dictionary can be updated in memory to run different simulations. An example script to simulate the model over the parameter ensemble is encoded in ``sample_ensemble.jl``. To execute this script, issue the command:

``julia> include("sample_ensemble.jl")``

This will load the parameter ensemble from disk (assumed to be stored in the ``data`` subdirectory), update the parameter dictionary returned by ``DataDictionary.jl`` with the new parameters, and solve the model equations. The Pareto rank of the simulation can be adjusted by changing the selection criteria (L13):

	# Select the desired rank -
	idx_rank = find(rank_array .<= 5.0)

This code will simulate all parameter sets in the ensemble with Pareto rank five or less.

The model ensemble was estimated using the [JuPOETs package](https://github.com/varnerlab/POETs.jl). The complement objective functions, and problem specific constraints are encoded in ``tf_network_lib.jl``. The JuPOETs package is described in the publication:

[Bassen D, Vilkhovoy M, Minot M, J. Butcher and J. Varner (2016) JuPOETs: A Constrained Multiobjective Optimization Approach to Estimate Biochemical Model Ensembles in the Julia Programming Language. bioRxiv 056044; doi: http://dx.doi.org/10.1101/056044](http://biorxiv.org/content/early/2016/05/30/056044)

__Prerequisites__: [Julia](http://julialang.org) and the [ODE package](https://github.com/JuliaDiffEq/ODE.jl) must be installed on your computer before the model equations can be solved. In addition, in the example routine ``sample_ensemble.jl`` the ensemble output is plotted using the [PyPlot](https://github.com/stevengj/PyPlot.jl) package which requires a working Python installation.  
