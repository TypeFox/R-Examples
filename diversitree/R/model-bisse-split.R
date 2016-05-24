## Split models should provide
##   1. make
##   2. make.cache
## Common other functions include:
##   branches.aux

## BiSSE/split is slightly more complicated than the other split
## models as dealing with the unresolved clades complicates matters a
## bit.

## 1: make
make.bisse.split <- function(tree, states, nodes, split.t=Inf,
                             unresolved=NULL, sampling.f=NULL,
                             nt.extra=10, strict=TRUE,
                             control=list()) {
  cache <- make.cache.bisse.split(tree, states, nodes, split.t,
                                  unresolved, sampling.f, nt.extra,
                                  strict)
  n.part <- cache$n.part
  unresolved <- cache$unresolved

  all.branches <- make.all.branches.split.dtlik(cache, control,
                                                initial.conditions.bisse)
  rootfunc <- make.rootfunc.split(cache, rootfunc.musse)

  ll <- function(pars, condition.surv=TRUE, root=ROOT.OBS,
                 root.p=NULL, intermediates=FALSE) {
    pars <- check.pars.bisse.split(pars, n.part)
    preset <- branches.unresolved.bisse.split(pars, unresolved)
    ans <- all.branches(pars, intermediates, preset)
    rootfunc(ans, pars, condition.surv, root, root.p, intermediates)
  }
  class(ll) <- c("bisse.split", "bisse", "dtlik", "function")
  ll
}

make.bisse.uneven <- function(tree, states, nodes, split.t=Inf,
                              unresolved=NULL, sampling.f=NULL,
                              nt.extra=10, strict=TRUE,
                              control=list()) {
  cache <- make.cache.bisse.split(tree, states, nodes, split.t,
                                  unresolved, sampling.f, nt.extra,
                                  strict)
  cache$info <- update.info.uneven(cache$info, make.info.bisse(tree))
    
  n.part <- cache$n.part
  unresolved <- cache$unresolved

  all.branches <- make.all.branches.split.dtlik(cache, control,
                                                initial.conditions.bisse)
  rootfunc <- make.rootfunc.split(cache, rootfunc.musse)

  ll <- function(pars, condition.surv=TRUE, root=ROOT.OBS,
                 root.p=NULL, intermediates=FALSE) {
    pars <- rep(list(check.pars.bisse(pars)), n.part)
    preset <- branches.unresolved.bisse.split(pars, unresolved)
    ans <- all.branches(pars, intermediates, preset)
    rootfunc(ans, pars, condition.surv, root, root.p, intermediates)
  }
  class(ll) <- c("bisse.uneven", "bisse", "dtlik", "function")
  ll
}

## Make requires the usual functions:
## 2: make.cache (initial.tip, root)
make.cache.bisse.split <- function(tree, states, nodes, split.t,
                                   unresolved, sampling.f,
                                   nt.extra, strict) {
  cache <- make.cache.bisse(tree, states, unresolved, NULL, nt.extra,
                            strict)
  cache <- make.cache.split.xxsse(tree, cache, nodes, split.t,
                                  sampling.f)
  unresolved <- cache$unresolved

  if ( !is.null(unresolved) ) {
    ## This ensures that the calculations should be slightly more
    ## identical by running out to the same number of species as for
    ## the non-split version.  Where one group mostly has smaller
    ## clades, this will slow the calculations down more than needed
    ## (may drop this after things settle down).
    n.part <- cache$n.part
    Nc.tot <- max(unresolved$Nc) + unresolved$nt.extra
    grp <- cache$group.branches[unresolved$target]
    tmp <- data.frame(unresolved[names(unresolved) != "nt.extra"],
                                   stringsAsFactors=FALSE)
    unresolved.split <- vector("list", n.part)
    
    for ( i in seq_len(n.part)[tabulate(grp, n.part) > 0] ) {
      j <- grp == i
      unresolved.split[[i]] <- as.list(tmp[j,])
      unresolved.split[[i]]$nt.extra <-
        Nc.tot - max(unresolved.split[[i]]$Nc)
    }

    cache$unresolved <- unresolved.split
  }

  cache
}

branches.unresolved.bisse.split <- function(pars, unresolved) {
  if ( is.null(unresolved) )
    return(NULL)
  ans <- list(target=integer(0), lq=numeric(0),
                 base=matrix(NA, 4, 0))
  for ( i in seq_along(pars) ) {
    unresolved.i <- unresolved[[i]]
    if ( !is.null(unresolved.i) ) {
      tmp <- branches.unresolved.bisse(pars[[i]], unresolved.i)
      ans$target <- c(ans$target, tmp$target)
      ans$lq <- c(ans$lq, tmp$lq)
      ans$base <- cbind(ans$base, tmp$base)
    }
  }
  ans
}

## Never used directly:
make.branches.aux.bisse <- function(cache, control)
  make.branches.aux.dtlik(cache, control)

check.pars.bisse.split <- function(pars, n.part) 
  check.pars.multipart(check.nonnegative(pars), n.part, 6)
