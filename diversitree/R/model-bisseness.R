## This file contains all additional functions necessary for running
## BiSSE-ness in R, as described in Magnuson-Ford, K., and
## S. Otto. 2012. Linking the investigations of character evolution
## and species diversification.  American Naturalist, XX: XX-XX.

## The BiSSE-ness model is identical to BiSSE in all aspects other
## than allowing character shifts at speciation.  This affects the
## following functions within diversitree.  The BiSSE analog of these
## functions have the same function name as those found here, omitting
## 'ness'.

## Models should provide:
##   1. make
##   2. info
##   3. make.cache, including initial tip conditions
##   4. initial.conditions(init, pars,t, idx)
##   5. rootfunc(res, pars, ...)

## Common other functions include:
##   stationary.freq
##   starting.point
##   branches

## TODO: Push this documentation around a bit.
## Using BiSSE-ness to estimate parameters (numbers follow other
## diversitree files)
##   1. make.bisseness is the overarching function that makes the
##      model
##
##   2. print.bisseness is the print 
##
##   3. argnames.bisseness gets and sets argument names
##
##   4. find.mle.bisseness alters the MLE search
##
##  (5. would be make.cache, but make.cache.bisse used)
##
##  (6. would be ll.bisseness, but done internally)
##
##   7. initial.conditions.bisseness combines the values from two
##      descending branches at a node, simultaneously accounting for
##      the possibility of cladogenetic change (eqn. 2 in the paper).
##
##   8. make.branches.bisseness computes calculation along a branch
##      using the differential equations specified in C file
##      'bisseness-eqs.c' (eqn. 3-6 in the paper).
##
##   9. branches.unresolved.bisseness calls the fortran code for
##      calculating the rate matrix and its exponent, to obtain the
##      probability of seeing the data listed in unresolved.
##
##  10. root.xxsseness specifies the calculations at the root, by
##      default conditioning the likelihood of the data on survival of
##      the two lineages descending from the root in a manner that
##      accounts for cladogenesis (eqn. 7; see also Nee et. al.,
##      1994).

## Using BiSSE-ness to simulate phylogenies:
##
##  11. stationary.freq.bisseness used for tree simulator to determine
##      the equilibrium root state; it can also be used to combine
##      probabilities at the root for parameter estimation using ML or
##      MCMC (but not the default, specify this explicitly).
##
##  12. tree.bisseness is the overall function used to simulate
##      phylogenies under BiSSE-ness.
##
##  13. make.tree.bisseness simulates phylogenies under BiSSE-ness.

## 1: make
make.bisseness <- function(tree, states, unresolved=NULL, sampling.f=NULL,
                           nt.extra=10, strict=TRUE, control=list()) {
  cache <- make.cache.bisseness(tree, states, unresolved, sampling.f,
                                nt.extra, strict)
  unresolved <- cache$unresolved
  if ( !is.null(cache$unresolved) )
    warning(paste("BiSSE-ness with unresolved clades has",
                  "not yet been extensively tested"))

  all.branches <- make.all.branches.dtlik(cache, control,
                                          initial.conditions.bisseness)
  rootfunc <- rootfunc.bisseness

  ll <- function(pars, condition.surv=TRUE, root=ROOT.OBS,
                 root.p=NULL, intermediates=FALSE) {
    check.pars.bisseness(pars)
    preset <- branches.unresolved.bisseness(pars, unresolved)
    ans <- all.branches(pars, intermediates, preset)
    rootfunc(ans, pars, condition.surv, root, root.p, intermediates)
  }
  class(ll) <- c("bisseness", "dtlik", "function")
  ll
}

## 2: info
make.info.bisseness <- function(phy) {
  list(name="bisseness",
       name.pretty="BiSSE-ness",
       ## Parameters:
       np=10L,
       argnames=default.argnames.bisseness(),
       ## Variables:
       ny=4L,
       k=2L,
       idx.e=1:2,
       idx.d=3:4,
       ## R version of derivative function:
       derivs=derivs.bisseness,
       ## Phylogeny:
       phy=phy,
       ## Inference:
       ml.default="subplex",
       mcmc.lowerzero=TRUE,
       ## These are optional
       doc=NULL,
       reference=c(
         "Magnuson-Ford and Otto (submitted)"))
}
default.argnames.bisseness <- function() 
  c("lambda0", "lambda1", "mu0", "mu1", "q01", "q10",
    "p0c", "p0a", "p1c", "p1a")

## 3: make.cache (mostly BiSSE's)
make.cache.bisseness <- function(tree, states, unresolved=NULL,
                                 sampling.f=NULL, nt.extra=10,
                                 strict=TRUE) {
  cache <- make.cache.bisse(tree, states, unresolved, sampling.f,
                            nt.extra, strict)
  cache$info <- make.info.bisseness(tree)
  cache
}

## 4: initial.conditions:
## Note that we ignore both 't' and 'idx'.
initial.conditions.bisseness <- function(init, pars, t, idx) {
  lambda <- pars[1:2]
  c(init[c(1,2),1],
    init[c(3,4),1] * init[c(3,4),2] * lambda * (1 - pars[c(7,9)]) +
    ##
    ((init[c(4,4),1] * init[c(3,3),2])/2 +
     (init[c(3,3),1] * init[c(4,4),2])/2) * lambda *
    c(pars[7]*pars[8], pars[9]*pars[10]) +
    ##
    init[c(4,3),1] * init[c(4,3),2] * lambda *
    c(pars[7]*(1-pars[8]), pars[9]*(1-pars[10])))
}

rootfunc.bisseness <- function(res, pars, condition.surv, root, root.p,
                               intermediates) {
  vals <- res$vals
  lq <- res$lq

  d.root <- vals[3:4]

  root.p <- root.p.calc(d.root, pars, root, root.p,
                        stationary.freq.bisseness)
  if ( condition.surv ) {
    lambda <- pars[1:2]
    e.root <- vals[1:2]

    nonextinct <-
      c(((1-pars[7])          * (1-e.root[1])^2+
         pars[7]*pars[8]      * (1-e.root[1])*(1-e.root[2]) +
         pars[7]*(1-pars[8])  * (1-e.root[2])^2),
        ((1-pars[9])          * (1-e.root[2])^2+
         pars[9]*pars[10]     * (1-e.root[1])*(1-e.root[2]) +
         pars[9]*(1-pars[10]) *(1-e.root[1])^2))
    
    d.root <- d.root / sum(root.p * lambda * nonextinct)
  }

  if ( root == ROOT.ALL )
    loglik <- log(d.root) + sum(lq)
  else
    loglik <- log(sum(root.p * d.root)) + sum(lq)

  if ( intermediates ) {
    res$root.p <- root.p
    attr(loglik, "intermediates") <- res
    attr(loglik, "vals") <- vals
  }

  loglik
}

branches.unresolved.bisseness <- function(pars, unresolved) {
  if ( is.null(unresolved) )
    return(NULL)
  Nc <- unresolved$Nc
  k <- unresolved$k
  nsc <- unresolved$nsc
  t <- unresolved$len
  nt <- max(Nc) + unresolved$nt.extra

  lambda0 <- pars[1]
  lambda1 <- pars[2]
  mu0 <- pars[3]
  mu1 <- pars[4]
  q01 <- pars[5]
  q10 <- pars[6]
  p0c <- pars[7]
  p0a <- pars[8]
  p1c <- pars[9]
  p1a <- pars[10]
  base <- nucexpl(nt, lambda0, lambda1, mu0, mu1, q01, q10,
                  p0c, p0a, p1c, p1a,
                  t, Nc, nsc, k)[,c(3,4,1,2),drop=FALSE]

  q <- rowSums(base[,3:4,drop=FALSE])
  base[,3:4] <- base[,3:4] / q

  ## Note the transpose here.
  list(target=unresolved$target,
       lq=log(q),
       base=t(base))
}

###########################################################################
## Additional functions
stationary.freq.bisseness <- function(pars) {
  lambda0 <- pars[1]
  lambda1 <- pars[2]
  mu0 <- pars[3]
  mu1 <- pars[4]
  q01 <- pars[5]
  q10 <- pars[6]
  p0c <- pars[7]
  p0a <- pars[8]
  p1c <- pars[9]
  p1a <- pars[10]

  g   <- (lambda0 - mu0) - (lambda1 - mu1)
  eps <- (lambda0 + mu0  +  lambda1 + mu1) * 1e-14
  ## ss=state shift: i.e. change within lineages *and* at speciation.
  ss0 <- q01 + lambda0*p0c*(2*(1 - p0a) + p0a)
  ss1 <- q10 + lambda1*p1c*(2*(1 - p1a) + p1a)
  
  if ( abs(g) < eps ) {
    if ( ss0 + ss1 == 0 ) 
      p <- 0.5
    else
      p <- ss1/(ss0 + ss1)
  } else {
    roots <- quadratic.roots(g, ss1 + ss0 - g, -ss1)
    roots <- roots[roots >= 0 & roots <= 1]
    if ( length(roots) > 1 )
      p <- NA
    else
      p <- roots
  }
  c(p, 1-p)
}

## For historical and debugging purposes, not used directly in the
## calculations, but branches function is generated this way
## internally.
make.branches.bisseness <- function(cache, control)
  make.branches.dtlik(cache$info, control)

check.pars.bisseness <- function(pars) {
  check.pars.nonnegative(pars, 10)
  if ( any(pars[7:10] > 1) ) # Already checked "< 0" above.
    stop("Probability parameters must lie between 0 and 1.")
  TRUE
}

derivs.bisseness <- function(t, y, pars) {
  E0 <- y[1]
  E1 <- y[2]
  D0 <- y[3]
  D1 <- y[4]

  la0 <- pars[1]
  la1 <- pars[2]
  mu0 <- pars[3]
  mu1 <- pars[4]
  q01 <- pars[5]
  q10 <- pars[6]
  p0c <- pars[7]
  p0a <- pars[8]
  p1c <- pars[9]
  p1a <- pars[10]

  c(-(mu0 + q01 + la0) * E0 + 
    la0 * E0 * E0 * (1 - p0c) + mu0 + q01 * E1 + 
    la0 * E0 * E1 * p0c * p0a + la0 * E1 * E1 * p0c * (1 - p0a),
    ##
    -(mu1 + q10 + la1) * E1 + 
    la1 * E1 * E1 * (1 - p1c) + mu1 + q10 * E0 + 
    la1 * E0 * E1 * p1c * p1a + la1 * E0 * E0 * p1c * (1 - p1a),
    ##
    -(mu0 + q01 + la0) * D0 + 
    2 * la0 * E0 * D0 * (1 - p0c) + q01 * D1 + 
    (E0 * D1 + E1 * D0) * la0 * p0c * p0a + 
    2 * la0 * E1 * D1 * p0c * (1 - p0a),
    ##
    -(mu1 + q10 + la1) * D1 + 
    2 * la1 * E1 * D1 * (1 - p1c) + q10 * D0 + 
    (E1 * D0 + E0 * D1) * la1 * p1c * p1a + 
    2 * la1 * E0 * D0 * p1c * (1 - p1a))
}
