# Multiple-reactivation model

### /!\ repository is under construction /!\

This repository contains jupyter notebooks, Stan models and data used in the paper

> [van Dorp *et al.*](https://insert.doi.here), Models of SIV rebound after treatment interruption that involve multiple reactivation events, submitted (2020)

To use the code, simply clone the repository

```bash
$ git clone https://github.com/lanl/multiple-reactivation-model.git
```

The model is described in full detail in the paper. In short, we present an refinement of a model developed by Pinkevych *et al.* and Hill *et al.* for SIV and HIV rebound after treatment interruption, that takes the effect of multiple reativation events into account. We derive a parameteric probability density function for the time to viral rebound, and implement this in the Stan programming language. We fit the model to data from SIV infected and very early treated macaques, published by Whitney *et al.* (2014), and (2018).

## Requirements

Make sure the following python packages are installed

- `pystan`
- `scipy`
- `numpy`
- `matplotlib`

For using the notebooks, Jupyter is required. 
Installation instructions can be found [here](https://jupyter.org/).
The code is only tested on Ubuntu 18.01, Python version 3.6.9 (pystan version 2.19), 
but *should* also work on Mac OSX and Windows 10.

## Data

The `data` folder contains the timeseries used in the analysis. The data has been described 
by Whitney *et al.* in two publications 

> [Whitney *et al.*](https://doi.org/10.1038/nature13594) Rapid seeding of the viral reservoir prior to SIV viraemia in rhesus monkeys. Nature, 512(7512):74-77, Aug 2014.

> [Whitney *et al.*](https://doi.org/10.1038/s41467-018-07881-9) Prevention of SIVmac251 reservoir seeding in rhesus monkeys by early antiretroviral therapy. Nat. Commun., 9(1):5429, 12 2018.

The data is separated into two csv files, containing rebound timeseries and acute infection timeseries respectively.

- `rebound-data-long.csv` Viral load measurements after treatment interruption
- `acute-data-long.csv` Viral load measurements during acute infection

## Notebooks

The following jupyter notebooks can be found in the `notebooks` folder:

- `VLProcessSimulations.ipynb` Simulate trajectories and first passage times of the viral load process. These are used to make **Figure 1 and 2** of the paper and **Figure S3, S4, S5 and S6** of the Supplementary Information.
- `ParseData.ipynb` This notebook is used to load the viral load data into a python dictionary. This dictionary is "pickled" into a file that can be loaded into the other notebooks for further analysis.
- ...
- `ReboundTimeVarianceAnalysis.ipynb` What percentage of the variance of the rebound time is due to the first recrudescence event, and when do other reactivations become important. This notebook is used to make **Figure 4**.
- ...

## Stan models

For more information on the Stan programming language, 
see the [Stan website](https://mc-stan.org/).
In the folder `stan-models` the following files can be found

- `logistic-rebound-model-multiple-reactiv-II.stan` is the main model, implementing the stochastic multiple-reactivation model in Stan. 
- `logistic-rebound-model-single-reactiv.stan` This is the simple single-reactivation model.
- `logistic-rebound-model-multiple-reactiv-CD.stan` This is [Pinkecych's approximation](https://dx.doi.org/10.1371%2Fjournal.ppat.1005740) of the multiple-reactivation model implemented in Stan.

## Python module

The Python module `mrmpytools` contains common functions used in the notebooks.
To get access to this module from a Jupyter notebook in the `notebooks` folder, 
add the parent directory to the path with

```py
import sys
sys.path.append("..") ## or the location of the mrmpytools module
```

Then import elements from this module. For instance 

```py
import mrmpytools.utilities as util
x1, x2 = util.unzip([('a', 1), ('b', 2)]) ## should be ['a', 'b'], [1, 2]
```

The folder `mrmpytools` contains the following files

- `plots.py` Some plotting functions
- `pystantools.py` Some functions for compiling and analyzing Stan models.
- ...

_______________________________________________________________


copyright 2020. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
Department of Energy/National Nuclear Security Administration. All rights in the program are
reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
Security Administration. The Government is granted for itself and others acting on its behalf a
nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare
derivative works, distribute copies to the public, perform publicly and display publicly, and to permit
others to do so.
