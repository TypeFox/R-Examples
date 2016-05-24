## Core ancestral state reconstruction (ASR) code.  I am implementing
## three different types of things:
##   (1) asr.marginal
##   (2) asr.joint
##   (3) asr.stoch
## There will perhaps be an overarching "asr" function some day but
## not right now.  All ASR functions take as a first argument a
## likelihood function, and dispatch based on the class of this.
## Methods currently implemented are mk2/mkn, but bisse and musse
## methods can be found in the unreleased package diversitree.unrel

## Core generics:
make.asr.marginal <- function(lik, ...) {
  UseMethod("make.asr.marginal")
}
make.asr.joint <- function(lik, ...) {
  UseMethod("make.asr.joint")  
}
make.asr.stoch <- function(lik, ...) {
  UseMethod("make.asr.stoch")
}

## Short cuts for one shot ASR:
asr.marginal <- function(lik, pars, nodes=NULL, ...)
  make.asr.marginal(lik)(pars, nodes, ...)
asr.joint <- function(lik, pars, n=1, ...)
  make.asr.joint(lik)(pars, n, ...)
asr.stoch <- function(lik, pars, n=1, ...)
  make.asr.stoch(lik)(pars, n, ...)

## Constrained functions require some care:
make.asr.marginal.constrained <- function(lik, ...) {
  lik.full <- attr(lik, "func")
  asr <- make.asr.marginal(lik.full, ...)
  function(pars, ...)
    asr(lik(pars, pars.only=TRUE), ...)
}

make.asr.joint.constrained <- function(lik, ...) {
  lik.full <- attr(lik, "func")
  asr <- make.asr.joint(lik.full, ...)
  function(pars, ...)
    asr(lik(pars, pars.only=TRUE), ...)
}

make.asr.stoch.constrained <- function(lik, ...) {
  lik.full <- attr(lik, "func")
  asr <- make.asr.stoch(lik.full, ...)
  function(pars, ...)
    asr(lik(pars, pars.only=TRUE), ...)
}

## Next, the utility functions for the different types of models This
## is to asr.marginal what all.branches is for the core models.  Here,
## the argument 'res' is the result of running all.branches
do.asr.marginal.R <- function(pars, cache, res, nodes, states.idx,
                              initial.conditions, branches, root,
                              ...) {
  ## Store these for easier calculation.
  children <- cache$children
  parent <- cache$parent
  len <- cache$len
  depth <- cache$depth
  root.idx <- cache$root
  anc <- cache$ancestors

  if ( is.null(nodes) )
    nodes <- root.idx:max(cache$order)
  else
    nodes <- nodes + cache$n.tip

  f <- function(nd) {
    ## Include current node but omit root:
    anc.nd <- c(nd, anc[[nd]])
    anc.nd <- anc.nd[-length(anc.nd)]
    p <- rep(NA, length(states.idx))

    for ( st in seq_along(states.idx) ) {
      lq <- res$lq
      branch.init <- res$init
      branch.base <- res$base
      branch.init[states.idx[-st],nd] <- 0
      y.in <- branch.init[,nd]

      for ( i in anc.nd ) {
        ans <- branches(y.in, len[i], pars, depth[i], i)
        lq[i] <- ans[[1]]
        branch.base[,i] <- ans[[2]]
        j <- parent[i]
        y.in <- initial.conditions(branch.base[,children[j,]], pars,
                                   depth[j], j)
        branch.init[,j] <- y.in
      }

      ans <- root(list(vals=branch.init[,root.idx], lq=lq), pars)

      ## explots IEEE arithmetic's exp(-Inf) == 0
      p[st] <- if ( is.na(ans) ) -Inf else ans
    }

    pp <- exp(p - max(p))
    pp/sum(pp)
  }

  matrix(unlist(lapply(nodes, f)), ncol=length(nodes))
}

make.do.asr.marginal <- function(all.branches, rootfunc) {
  eb <- environment(all.branches)
  cache <- eb$cache
  states.idx <- cache$info$idx.d
  if (isTRUE(cache$info$partitioned)) {
    branches <- eb$branches.split
    initial.conditions <- eb$initial.conditions.split
  } else {
    branches <- eb$branches
    initial.conditions <- eb$initial.conditions
  }
  function(pars, nodes, preset, ...) {
    root.f <- function(res, pars)
      rootfunc(res, pars, ...)
    res <- all.branches(pars, TRUE, preset)
    do.asr.marginal.R(pars, cache, res, nodes, states.idx,
                      initial.conditions, branches, root.f)
  }
}


## Utility function for drawing one or more samples from the joint
## distribution.

## li is k * len matrix; for a node n, li[,n] comes in the order
##   Pr(D_n|1), Pr(D_n|2), ..., Pr(D_n|k)
## Pr(D|i) is the conditional probability of the data conditional on a
## node being in the state 'i'.  It .

## pij is a (k*k) * len column matrix; the column n comes in the order
##   p11, p21, ..., pk1, p12, ..., pkk
## so that
##   matrix(pij[,nd], k, k)
## is a matrix with where m[i,j] is the probability of moving from
## state i to state j.
do.asr.joint <- function(n, k, order.C, parent.C, li, pij, root.p,
                         as.01) {
  if ( n > 1 ) {
    ret <- matrix(NA, n, length(order.C))
    for ( i in seq_len(n) )
      ret[i,] <- .Call("r_do_asr_joint", k, order.C, parent.C, li,
                       pij, root.p, as.01, PACKAGE="diversitree")
    ret
  } else {
    .Call("r_do_asr_joint", k, order.C, parent.C, li,
          pij, root.p, as.01, PACKAGE="diversitree")
  }
}

