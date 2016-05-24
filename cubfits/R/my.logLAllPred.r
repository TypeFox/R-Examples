### This either returns a logL vector or sum of the vector.
###
### These functions are for all predicting genes.

### Get the specific function according to the options.
get.my.logLAllPred <- function(model.Phi){
  if(!any(model.Phi[1] %in% .CF.CT$model.Phi)){
    stop("model.Phi is not found.")
  }
  ret <- eval(parse(text = paste("my.logLAllPred.",
                                 model.Phi[1], sep = "")))
  assign("my.logLAllPred", ret, envir = .cubfitsEnv)
  ret
} # End of get.my.logLAllPred().


### Function to calculate complete logL for
### (phi, b) given y and n
my.logLAllPred.lognormal <- function(phi, y, n, b, reu13.df = NULL){
  .cubfitsEnv$my.logdmultinomCodAllR(b, phi, y, n, reu13.df = reu13.df)
} # End of my.logLAllPred.lognormal().

### No need to changed from my.logLAll.lognormal since prior does not count.
my.logLAllPred.logmixture <- function(phi, y, n, b,
    reu13.df = NULL){
  my.logLAllPred.lognormal(phi, y, n, b, reu13.df = reu13.df)
} # End of my.logLAll.logmixture().

