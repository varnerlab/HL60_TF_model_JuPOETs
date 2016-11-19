## HL60 Gene Regulatory Network Model


### Background ###
In this study, we present an effective model All-Trans Retinoic Acid (ATRA)-induced differentiation of HL-60 cells.
The model describes a key architectural feature of ATRA-induced differentiation, positive feedback between an ATRA-inducible signalsome complex, the activation of the mitogen activated protein kinase (MAPK) cascade and the downstream gene expression program.
The model, which was developed by integrating logical rules with kinetic modeling, was significantly smaller than previous HL60 models.
However, despite its simplicity, it captured key features of ATRA induced differentiation of HL-60 cells.
We identified an ensemble of effective model parameters using measurements taken from ATRA-induced HL-60 cells.
Using these parameters, we were able to simulate key features of the signaling and gene expression programs following from ATRA addition.

### Installation

You can download this repository as a zip file, or clone or pull it by using the command:

	git pull https://github.com/varnerlab/HL60_TF_model_JuPOETs.git

or

	git clone https://github.com/varnerlab/HL60_TF_model_JuPOETs.git

### Model code and parameter ensemble
The HL60 GRN equations were implemented in [Julia](http://julialang.org) and solved using the ode23s routine of the [ODE package](https://github.com/JuliaDiffEq/ODE.jl). The model code and parameter ensemble is freely available under an [MIT software license](https://opensource.org/licenses/MIT).

The model equations are encoded in ``Balances.jl`` which is called by the ``SolveBalances.jl`` driver function. The user should not directly call ``SolveBalances.jl``. Rather, multiple parameter sets can be simulated by calling the driver function from a script. The kinetic and other model parameters are encoded in ``DataDictionary.jl`` as a dictionary. The parameters stored in this dictionary can be updated in memory to run different simulations. An example script to simulate the model over the parameter ensemble is encoded in ``sample_best_ensemble.jl``. To execute this script, issue the command:

``julia> include("sample_best_ensemble.jl")``

This will load the parameter ensemble from disk, update the parameter dictionary returned by ``DataDictionary.jl`` with the new parameters, and solve the model equations. The model solution is written to the ``best_ensemble`` directory.

The model ensemble was estimated using the ``ReallySimpleSearch.jl`` routine. For details about the model, the order of the species or the order of the parameters, please look at ``DataDictionary.jl``.

__Prerequisites__: [Julia](http://julialang.org) and the [ODE package](https://github.com/JuliaDiffEq/ODE.jl) must be installed on your computer before the model equations can be solved. In addition, in the example routine ``sample_best_ensemble.jl`` the ensemble output is plotted using the [PyPlot](https://github.com/stevengj/PyPlot.jl) package which requires a working Python installation.  
