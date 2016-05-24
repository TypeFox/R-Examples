################################################
#### AUTHOR:     Arnost Komarek             ####
####             (2005)                     ####
####                                        ####
#### FILE:       bayessurvreg2.priorInit.R  ####
####                                        ####
#### FUNCTIONS:  bayessurvreg2.priorInit    ####
################################################

### ======================================
### bayessurvreg2.priorInit
### ======================================
## Subsection of bayessurvreg2
## -> just to make it more readable
##
## 30/01/2005
##
## REMARK: both design$Y and design2$Y contain not-transformed times (i.e. t and not y)
##
bayessurvreg2.priorInit <- function(prior, init, design, mcmc.par, prior2, init2, design2, mcmc.par2, doubly)
{
  return(bayesBisurvreg.priorInit(dim=1, prior=prior, init=init, design=design, mcmc.par=mcmc.par,
                                         prior2=prior2, init2=init2, design2=design2, mcmc.par2=mcmc.par2, doubly=doubly))
}
