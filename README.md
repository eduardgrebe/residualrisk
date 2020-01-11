# Weusten et al. (2011) model for residual risk of HIV transfusion transmission with NAT screening

This code implements a mathematical model for deriving "risk day equivalents" (i.e. infectious window period), in the case of nucleic acid screening of blood donations for viruses like HIV, for the purposes of estimating residual risk of transfusion transmission. In this implementation, we assume that the user will have an estimate of incidence in the donor population, which when multiplied by the "risk day equivalents" estimate will yield the risk that any given donation will transmit the virus.

The model is described in:
* Weusten, J., Vermeulen, M., van Drimmelen, H. and Lelie, N. (2011), Refinement of a viral transmission risk model for blood donations in seroconversion window phase screened by nucleic acid testing in different pool sizes and repeat test algorithms. Transfusion, 51: 203-215. doi:[10.1111/j.1537-2995.2010.02804.x](https://doi.org/10.1111/j.1537-2995.2010.02804.x).

Code copyright Eduard Grebe <EGrebe@vitalant.org>, 2019. Released under the The GNU Affero General Public License Version 3.0.

## Suggested citation

If you use this code in published research, I suggest citing both the Weusten article and this repository, by means of the DOI on Zenodo:
* Grebe, E. (2019) Residual risk of HIV transfusion transmission with NAT screening. Jupyter Notebook. Version 1.0. doi:[10.5281/zenodo.3588570](https://doi.org/10.5281/zenodo.3588570).


## Running

The code is written for Python 3.

The following non-standard python modules are required:
* numpy
* scipy
* jupyterlab (to render the notebook; the classic Jupyter notebook server should also work)
* matplotlib (to render the plots in the example)

If you use conda, you can set up an environment called `residualrisk` by using the environment specification included with this repository and activate it by executing:

```
conda env create -f environment.yml
conda activate residualrisk
```

Or you can install the dependencies with pip.
