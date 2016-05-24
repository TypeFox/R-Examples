################################################################################
#  Baseline hazard distributions                                               #
################################################################################
#                                                                              #
#  These are functions with parameters                                         #
#   - pars: the vector of parameters                                           #
#   - t   : the time point                                                     #
#   - what: the quantity to be returned by the function                        #
#           either  "H", the cumulated hazard,                                 #
#               or "lh", the log-hazard                                        #
#                                                                              #
#                                                                              #
#   Date:                 December, 19, 2011                                   #
#   Last modification on: June, 27, 2012                                       #
################################################################################


################################################################################
################################################################################


################################################################################
#                                                                              #
#   Weibull baseline hazard function                                           #
#                                                                              #
#   Parameters:                                                                #
#    [1] rho    > 0                                                            #
#    [2] lambda > 0                                                            #
#                                                                              #
#   Hazard:                                                                    #
#    h(t) = \rho \lambda t^(\rho-1)                                            #
#                                                                              #
#                                                                              #
#   Date: December, 19, 2011                                                   #
################################################################################

weibull <- function(pars,
                    t, 
                    what){
  if (what == "H")
    return(pars[2] * t^(pars[1]))
  else if (what == "lh")
    return(log(pars[1]) + log(pars[2]) + ((pars[1] - 1) * log(t)))
}


################################################################################
#                                                                              #
#   Inverse Weibull baseline hazard function                                   #
#                                                                              #
#   Parameters:                                                                #
#    [1] rho    > 0                                                            #
#    [2] lambda > 0                                                            #
#                                                                              #
#   Hazard:                                                                    #
#    h(t) = [\rho / (\lambda t^(\rho + 1))] / [exp{1 /(\lambda t^\rho)} - 1]   #
#    H(t) = log[exp{1 /(\lambda t^\rho)} - 1]                                  #
#                                                                              #
#                                                                              #
#   Date:                 June, 26, 2012                                       #
#   Last modification on: June, 27, 2012                                       #
################################################################################

inweibull <- function(pars,
                       t, 
                       what){
  if (what == "H")
    return(-log(1 - exp(-1 / (pars[2] * t^(pars[1])))))
  else if (what == "lh")
    return(log(pars[1]) - log(pars[2]) - ((pars[1] + 1) * log(t)) -
      log(exp(1 / (pars[2] * t^(pars[1]))) - 1))
}


################################################################################
#                                                                              #
#   Exponential baseline hazard function                                       #
#                                                                              #
#   Parameters:                                                                #
#    [1] lambda > 0                                                            #
#                                                                              #
#   Hazard:                                                                    #
#    h(t) = \lambda                                                            #
#                                                                              #
#                                                                              #
#   Date: December, 19, 2011                                                   #
################################################################################

exponential <- function(pars,
                        t, 
                        what){
  if (what == "H")
    return(pars * t)
  else if (what == "lh") 
    return(log(pars))
}



################################################################################
#                                                                              #
#   Gompertz baseline hazard function                                          #
#                                                                              #
#   Parameters:                                                                #
#    [1] gamma  > 0                                                            #
#    [2] lambda > 0                                                            #
#                                                                              #
#   Hazard:                                                                    #
#    h(t) = \lambda exp(\gamma t)                                              #
#                                                                              #
#                                                                              #
#   Date: December, 19, 2011                                                   #
################################################################################

gompertz <- function(pars,
                     t, 
                     what){
  if (what == "H") 
    return(pars[2] / pars[1] * (exp(pars[1] * t) - 1))
  else if (what == "lh") 
    return(log(pars[2]) + pars[1] * t)
}



################################################################################
#                                                                              #
#   Lognormal baseline hazard function                                         #
#                                                                              #
#   Parameters:                                                                #
#    [1] mu     \in \mathbb R                                                  #
#    [2] sigma  > 0                                                            #
#                                                                              #
#   Hazard:                                                                    #
#    h(t) = \frac{ \phi(\log t -\mu \over \sigma) }{                           #
#             \sigma t [1-\Phi(\log t -\mu \over \sigma)] }                    #
#                                                                              #
#   with \phi the density and \Phi the distribution function                   #
#   of a standard Normal                                                       #
#                                                                              #
#                                                                              #
#   Date: December, 19, 2011                                                   #
################################################################################

lognormal <- function(pars,
                      t, 
                      what){
  if (what == "H")  #return - log (S)
    return(- log(1 - plnorm(t, meanlog=pars[1], sdlog=pars[2])))
  else if (what == "lh")  #return log(f) - log(S)
    return(dlnorm(t, meanlog=pars[1], sdlog=pars[2], log=TRUE) -
      log(1 - plnorm(t, meanlog=pars[1], sdlog=pars[2])))
}



################################################################################
#                                                                              #
#   Loglogistic baseline hazard function                                       #
#                                                                              #
#   Parameters:                                                                #
#    [1] alpha \in \mathbb R                                                   #
#    [2] kappa > 0                                                             #
#                                                                              #
#   Hazard:                                                                    #
#    h(t) = \frac{ \exp(\alpha) \kappa t^{\kappa-1} }{                         #
#             1 + \exp(\alpha) t^\kappa }                                      #
#                                                                              #
#                                                                              #
#                                                                              #
#   Date: December, 19, 2011                                                   #
################################################################################

loglogistic <- function(pars,
                        t, 
                        what){
  if (what == "H") 
    return(log(1 + exp(pars[1]) * t^(pars[2])))
  else if (what == "lh") 
    return(pars[1] + log(pars[2]) + (pars[2] - 1) * log(t) - 
      log(1 + exp(pars[1]) * t^pars[2]) )
}
