###################################################
#### AUTHOR:     Arnost Komarek                ####
####             (2005)                        ####
####                                           ####
#### FILE:       bayessurvreg3.priorBeta.R     ####
####                                           ####
#### FUNCTIONS:  bayessurvreg3.priorBeta       ####
###################################################

### ======================================
### bayessurvreg3.priorBeta
### ======================================
bayessurvreg3.priorBeta <- function(prior.beta, init, design)
{
  return(bayessurvreg2.priorBeta(prior.beta, init, design))
}
