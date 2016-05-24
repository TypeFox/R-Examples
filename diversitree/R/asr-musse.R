## (1) Marginal ASR:
## This should work without modification for:
##   musse
##   musse.split
##   musse.t
##   musse.td
make.asr.marginal.musse <- function(lik, ...) {
  e <- environment(lik)
  f.pars <- e$f.pars
  do.asr <- make.do.asr.marginal(e$all.branches, e$rootfunc)
  asr <- function(pars, nodes=NULL, condition.surv=TRUE,
                  root=ROOT.FLAT, root.p=NULL, ...) {
    pars2 <- f.pars(pars)
    do.asr(pars2, nodes, NULL, # below here extra args to rootfunc:
           condition.surv, root, root.p, intermediates=FALSE)
  }
  asr
}
