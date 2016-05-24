##generic
DIC <- function(mod, return.pD = FALSE, ...) {
  UseMethod("DIC", mod)
}



##default
DIC.default <- function(mod, return.pD = FALSE, ...) {
  stop("\nFunction not yet defined for this object class\n")
}


##bugs
DIC.bugs <- function(mod, return.pD = FALSE, ...){
  ##  DIC = posterior mean of the deviance + pD, where pD is the effective number of parameters
  if(return.pD == FALSE){
    DIC <- mod$DIC
  } else {DIC <- mod$pD}
  return(DIC)
}



##jags
DIC.rjags <- function(mod, return.pD = FALSE, ...){
  ##  DIC = posterior mean of the deviance + pD, where pD is the effective number of parameters
  if(return.pD == FALSE){
    DIC <- mod$BUGSoutput$DIC
  } else {DIC <- mod$BUGSoutput$pD}
  return(DIC)
}
