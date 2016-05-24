# computing DCC

dcc.est <- function(dvar, param){
   uncR <- cov(dvar)
   out <- .Call("dcc_est", dvar, uncR, param[1], param[2])
   list(DCC=out[[1]], Q=out[[2]])
}