## About

This application implements a model for the infectious window period (IWP) with 
NAT screening for a viral transfusion-transmissible infection. The default 
parameters are appropriate for HIV, but with appropriate parameters the logic is
equally applicable to other infections like HBV and HCV. The infectious-but-
undetectable window period depends on:
* the probability of transfusion transmission as a function of viral 
concentration in the transfused component
* the probability of non-detection as a function of viral concentration in the
screening sample (which depends on the sensitivity of the NAT screening assay,
minipool size and the retesting procedure)
* the viral growth rate (or doubling time) during early infection ("ramp-up
viremia")

The infectious window period is given by the area under the joint probability of
transfusion-transmission and non-detection.

The residual risk of transfusion transmission despite screening is given by the
probability that a random donation will be infection-positive and very recently
acquired -- i.e. incidence x infectious window period.

The model incorporates a "dose-response" model for the probability of 
transfusion-transmission:
* Belov, AA et al. (2020) Transfusion. Submitted.

and an adaptation of the Weusten model for infectious window period:
* Weusten, J. et al. (2011) Transfusion.

Written by Eduard Grebe <Eduard.Grebe@ucsf.edu>.

This program is free software: you can redistribute it and/or modify
it under the terms of the [GNU Affero General Public License]
(https://www.gnu.org/licenses/agpl-3.0.en.html) as published by the Free 
Software Foundation, either version 3 of the License, or (at your option) any 
later version.
