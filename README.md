# JEVFellowshipAnnex

This repository contains the files needed to run the analyses in the research proposal: 'From alternative hosts to alternative interventions: identifying drivers of epidemic dynamics for Japanese encephalitis virus in Bangladesh', submitted to the MRC Career Development Award on 21st April 2021.

The code is made available under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>. You are free to reuse this code provided that you give appropriate credit, provide a link to the license, and indicate if changes were made. You may do so in any reasonable manner, but not in any way that suggests the licensor endorses you or your use. 

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />


## Sample size estimates for bird-baited traps
- **r_mozBirdSampling.R** - this R script contains the simulations used to estimate the number of replicates required to determine differences in the number of mosquitoes collected from three different bird species in two different locations.

## Sample size estimates for mosquito blood feeding
- **r_mozBloodmealSampling.R** - this R script contains the simulations used to determine the sampling strategy to quantify the effect of household host community composition on the proportion of vector blood meals taken on each host species present.

## Sample size estimates for estimates of JEV prevalence in mosquitoes
- **r_odeThreeHostVarMoz.R** - this file contains an ordinary differential equation model of Japanese encephalitis virus transmission used to help inform the magnitude of possible seasonal changes in mosquito JEV prevalence.

- **r_odeThreeHostVarMozDatSimSensitivity.R** - this file contains the code requried to carry out a sensitivity analysis of the ODE model.

- **r_odeThreeHostVarMozSummaryFuncs.R** - this file contains summary functions used for plotting ODE model outputs

- **r_MLEMozPrevFuncs.R** - this file contains the functions to estimate the maximum likelihood prevalence of JEV in mosquitoes from pooled samples from data. It also contains the code required to simulate data.

- **r_MLEMozPrevalencePwr.R** - this file contains the simulations required to estimate sample size for estimating mosquito JEV prevalence.

## Sample size estimates for livestock serosurveys
- **r_MLEforceOfInfectionFuncs.R** - this file contains the code required to estimate the force of infection from mosquitoes to livestock based on age-structured seroprevalence data

- **r_MLEforceofInfectionPwr.R** - this file contains the simulations to estimate sample size for livestock force of infection.

- **r_livestockSeroprevalenceSampleSize.R** - this file contains the simulations to estimate sample size for livestock seroprevalence.
