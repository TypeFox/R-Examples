make.asr.marginal.geosse <- function(lik, ...) {
  e <- environment(lik)
  do.asr <- make.do.asr.marginal(e$all.branches, e$rootfunc)
  asr <- function(pars, nodes=NULL, condition.surv=TRUE,
                  root=ROOT.FLAT, root.p=NULL, ...) {
    check.pars.nonnegative(pars, 7)    
    do.asr(pars, nodes, NULL, # below here extra args to rootfunc:
           condition.surv, root, root.p, intermediates=FALSE)
  }
  asr
}
