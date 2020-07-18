# Residual risk of HIV transfusion transmission with NAT screening

*Eduard Grebe <EGrebe@vitalant.org>, 17 July 2020*

The code in this repository implements a tool for estimating the residual risk 
of HIV transfusion transmission despite screening of blood donation samples with 
nucleic acid testing (NAT). The central concepts underlying the approach is the 
"infectious window period" (IWP) and HIV incidence in the donor population -- 
reflecting the probability that a random donation will be in the infectious 
window period. The infectious window period is a function of the viral growth 
rate during early ramp-up viraemia, the sensitivity of the screening assay and 
the volume transfused. Multiplying the infectious window period by incidence 
yields the probability (per transfusion) of a transfusion transmission.

This "incidence times window period" model is well-established and is described 
in:

* Busch MP, Glynn SA, Stramer SL, Strong DM, Caglioti S, Wright DJ, Pappalardo B, Kleinman SH. A new strategy for estimating risks of transfusion-transmitted viral infections based on rates of detection of recently infected donors. Transfusion. 2005;45:254–264. doi:[10.1111/j.1537-2995.2004.04215.x](https://doi.org/10.1111/j.1537-2995.2004.04215.x).

Busch et al. presumed infectivity at a concentration of 1 RNA copy/20mL (i.e. 
the approximate volume of packed red blood cell components). They estimated the 
window period from presumed infectivity to individual donation NAT (ID-NAT) 
detection at 5.6 days, and from ID-NAT to minipool NAT (MP-NAT) with minipools 
of 16 samples at 3.4 days.

Their approach was an application of the general logic described Weusten and 
colleagues in in 2002:

* Weusten JJAM, Van Drimmelen HAJ, Lelie PN. Mathematic modeling of the risk of HBV, HCV, and HIV transmission by window-phase donations not detected by NAT. Transfusion. 2002;42:537–548.

The same first author led a refinement of this mathematic model published in 
2011, and on which this implementation draws heavily:

* Weusten JJAM, Vermeulen M, Van Drimmelen HAJ, Lelie PN. Refinement of a viral transmission risk model for blood donations in seroconversion window phase screened by nucleic acid testing in different pool sizes and repeat test algorithms. Transfusion. 2011;51:203-215. doi:[10.1111/j.1537-2995.2010.02804.x](https://doi.org/10.1111/j.1537-2995.2010.02804.x).

This software tool implementes the same basic logic as the Weusten et al. (2011) 
model, but diverges from it in the following important ways:

1. the probability of a transfusion transmission as a function of time since infection is modeled using a log-linear viral growth model (Fiebig et al., 2003) and an HIV-TT dose response model by Belov et al. (submitted), that was fit to data from a non-human primate SIV transmission study (Ma et al., 2009);
2. we are agnostic as to the source of HIV incidence estimates, which are treated as external to the model, rather than relying on repeat donor seroconversion patterns; and
3. we provide a method for conducting an uncertainty analysis by specifying distributions around critical model parameters and drawing sets of parameter combinations to produce credible ranges around infectious window period and residual risk estimates.

It also draws inspiration from a residual risk estimation tool previously 
developed and released by Alex Welte, Eduard Grebe and Andrew Powrie as part of 
an HIV infection dating tool:

* Grebe E, Facente SN, Powrie A, Gerber J, Grootboom K, Priede G, Welte A. Infection Dating Tool. Zenodo. 2018. doi:[10.5281/zenodo.1488117](https://doi.org/10.5281/zenodo.1488117).

This tool is a work and progress and contributions are encouraged. Enhancements 
currently contemplated include incorporating enhancements to the underlying 
model by Jeremy Bingham under the supervision of Alex Welte.

The code in the Jupyter Notebook and R-Shiny web application was written by 
Eduard Grebe, and the dose response model fitting code in the RMarkdown notebook 
was written by Artur Belov.

Copyright is asserted on the computer code, and it is released under the terms 
of the GNU Affero General Public License Version 3.0. You are free to reproduce, 
modify and distribute copies of it under the terms of this license.

## Suggested citation

If you use this code in published research, we suggest citing both the Weusten 
article and this repository, by means of the DOI on Zenodo:

* Grebe, E. (2019) Residual risk of HIV transfusion transmission with NAT screening. Version 1.0. doi:[10.5281/zenodo.3588570](https://doi.org/10.5281/zenodo.3588570).

**Version 1.5 is in progress.**

## Running

The core model code is written in Python 3. The web application is written in R 
using the Shiny framework. The R code for the dose response model requires 

The core Python code requires the following modules:
* numpy
* scipy

To edit and excute the Jupyter notebook requires:

* jupyterlab (to render the notebook; the classic Jupyter notebook server should also work)
* matplotlib (to render the plots in the example)

If you use conda, you can set up an environment called `residualrisk` by using the environment specification included with this repository and activate it by executing:

```
conda env create -f environment.yml
conda activate residualrisk
```

Or you can install the dependencies with pip.
