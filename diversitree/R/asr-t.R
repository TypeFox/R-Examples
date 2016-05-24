make.asr.marginal.dtlik.t <- function(lik, ...) {
  e <- environment(lik)
  cache <- get.cache(lik)
  do.asr <- make.do.asr.marginal(e$all.branches, e$rootfunc)
  preset <- NULL

  ## Slightly different prototypes will be needed for different
  ## models, I guess.  A mkn model doesn't want condition.surv.  In
  ## contrast with other models, no parameter checking here (happens
  ## automatically)
  asr <- function(pars, nodes=NULL, condition.surv=TRUE,
                  root=ROOT.FLAT, root.p=NULL)
    do.asr(pars, nodes, preset, # below here extra args to root func
           condition.surv, root, root.p, intermediates=FALSE)

  asr
}

make.asr.marginal.bisse.t <- make.asr.marginal.dtlik.t
make.asr.marginal.musse.t <- make.asr.marginal.dtlik.t

