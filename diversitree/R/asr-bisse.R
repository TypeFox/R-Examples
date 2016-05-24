## (1) Marginal ASR:
## This should work without modification for:
##   bisse
##   bisse.split
##   bisse.t
##   bisse.td
make.asr.marginal.bisse <- function(lik, ...) {
  e <- environment(lik)
  unresolved <- e$unresolved
  do.asr <- make.do.asr.marginal(e$all.branches, e$rootfunc)
  asr <- function(pars, nodes=NULL, condition.surv=TRUE,
                  root=ROOT.FLAT, root.p=NULL, ...) {
    check.pars.bisse(pars)
    preset <- branches.unresolved.bisse(pars, unresolved)
    do.asr(pars, nodes, preset, # below here extra args to rootfunc:
           condition.surv, root, root.p, intermediates=FALSE)
  }
  asr
}

make.asr.marginal.bisse.split <- function(lik, ...) {
  e <- environment(lik)
  unresolved <- e$unresolved
  cache <- get.cache(lik)
  do.asr <- make.do.asr.marginal(e$all.branches, e$rootfunc)
  asr <- function(pars, nodes=NULL, condition.surv=TRUE,
                  root=ROOT.FLAT, root.p=NULL, ...) {
    pars <- check.pars.bisse.split(pars, cache$n.part)
    preset <- branches.unresolved.bisse.split(pars, unresolved)
    do.asr(pars, nodes, preset, # below here extra args to rootfunc:
           condition.surv, root, root.p, intermediates=FALSE)
  }
  asr
}
