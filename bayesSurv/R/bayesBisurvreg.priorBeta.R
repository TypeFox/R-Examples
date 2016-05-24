#################################################
#### AUTHOR:     Arnost Komarek              ####
####             (2005)                      ####
####                                         ####
#### FILE:       bayesBisurvreg.priorBeta.R  ####
####                                         ####
#### FUNCTIONS:  bayesBisurvreg.priorBeta    ####
#################################################

### ======================================
### bayesBisurvreg.priorBeta
### ======================================
## Manipulation with the prior specification for the regression parameters
## (and means of random effects)
## - version for bayesBisurvreg
##
## 14/01/2005
## 24/01/2005: changed to a call to function bayessurvreg2.priorBeta
## =========================================================================
bayesBisurvreg.priorBeta <- function(prior.beta, init, design)
{
  return(bayessurvreg2.priorBeta(prior.beta, init, design))
} 
