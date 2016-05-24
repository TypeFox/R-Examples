## GeoSSE model, by Emma Goldberg <eeg@uic.edu>

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

## 1: make
make.geosse <- function(tree, states, sampling.f=NULL, strict=TRUE,
                        control=list()) {
  cache <- make.cache.geosse(tree, states, sampling.f, strict)

  all.branches <- make.all.branches.dtlik(cache, control,
                                          initial.conditions.geosse)
  rootfunc <- rootfunc.geosse

  ll <- function(pars, condition.surv=TRUE, root=ROOT.OBS,
                 root.p=NULL, intermediates=FALSE) {
    check.pars.geosse(pars)
    ans <- all.branches(pars, intermediates)
    rootfunc(ans, pars, condition.surv, root, root.p, intermediates)
  }
  class(ll) <- c("geosse", "dtlik", "function")
  ll
}

## 2: info
make.info.geosse <- function(phy) {
  list(name="geosse",
       name.pretty="GeoSSE",
       ## Parameters:
       np=7L,
       argnames=default.argnames.geosse(),
       ## Variables:
       ny=6L,
       k=3L,
       idx.e=1:3,
       idx.d=4:6,
       ## R version of derivative function:
       derivs=derivs.geosse,
       ## Phylogeny:
       phy=phy,
       ## Inference:
       ml.default="subplex",
       mcmc.lowerzero=TRUE,
       ## These are optional
       doc=NULL,
       reference=c(
         "Goldberg et al. (2011) doi:10.1093/sysbio/syr046"))
}
default.argnames.geosse <- function()
  c("sA", "sB", "sAB", "xA", "xB", "dA", "dB")

## 3: make.cache (& initial.tip)
## Note: classe uses the same functions as musse here to generate the
## cache and initial tip states, though we have to move from 0..2 to
## 1..3 state indexing to make this work smoothly.
make.cache.geosse <- function(tree, states, sampling.f=NULL,
                              strict=TRUE) {
  cache <- make.cache.musse(tree, states+1L, 3L, sampling.f, strict)
  cache$info <- make.info.geosse(tree)
  cache
}

## 4: initial.conditions:
## Note that we ignore both 't' and 'idx'.
initial.conditions.geosse <- function(init, pars, t, idx) {
  ## E.0, E.1, E.2
  e <- init[c(1,2,3),1]

  ## D.1, D.2  (Eq. 6bc)
  d12 <- init[c(5,6),1] * init[c(5,6),2] * pars[c(1,2)]

  ## D.0 (Eq. 6a)
  d0 <- 0.5 * sum(init[c(4,5),1] * init[c(5,4),2] * pars[1] + 
                  init[c(4,6),1] * init[c(6,4),2] * pars[2] +
                  init[c(5,6),1] * init[c(6,5),2] * pars[3])
  d <- c(d0, d12)

  c(e, d)
}

## 5: rootfunc:
rootfunc.geosse <- function(res, pars, condition.surv, root, root.p,
                            intermediates) {
  vals <- res$vals
  lq <- res$lq

  d.root <- vals[4:6]
  root.p <- root.p.calc(d.root, pars, root, root.p,
                        stationary.freq.geosse)
  if ( condition.surv ) {
    e.root <- vals[1:3]
    ## AB species are subject to all three speciation rates  
    lambda <- c(sum(pars[1:3]), pars[1:2])
    d.root <- d.root/sum(root.p * lambda * (1 - e.root)^2)
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

###########################################################################
## Additional functions
stationary.freq.geosse <- function(pars) {
  sA  <- pars[1]
  sB  <- pars[2]
  sAB <- pars[3]
  xA  <- pars[4]
  xB  <- pars[5]
  dA  <- pars[6]
  dB  <- pars[7]

  A <- matrix(c(
               -xA - xB - sAB,   dA,             dB,
                sA + xB + sAB,   sA - dA - xA,   0,
                sB + xA + sAB,   0,              sB - dB - xB
               ), nrow=3, byrow=TRUE)

  ## continuous time, so the dominant eigenvalue is the largest one
  ## return its eigenvector, normalized  
  evA <- eigen(A)
  i <- which(evA$values == max(evA$values))
  evA$vectors[,i] / sum(evA$vectors[,i])
}

## Replacement function from Emma, 4 Nov 2010.
## from Magallon & Sanderson (2001), rather than a bd fit
starting.point.geosse <- function(tree, eps=0.5) {
  if (eps == 0) {
    s <- (log(Ntip(tree)) - log(2)) / max(branching.times(tree))
    x <- 0
    d <- s/10
  } else {
    n <- Ntip(tree)
    r <- ( log( (n/2) * (1 - eps*eps) + 2*eps + (1 - eps)/2 *
               sqrt( n * (n*eps*eps - 8*eps + 2*n*eps + n))) - log(2)
          ) / max(branching.times(tree))
    s <- r / (1 - eps)
    x <- s * eps
    d <- x
  }
  p <- c(s, s, s, x, x, d, d)
  names(p) <- default.argnames.geosse()
  p
}

check.pars.geosse <- function(pars)
  check.pars.nonnegative(pars, 7)

## For historical and debugging purposes, not used directly in the
## calculations, but branches function is generated this way
## internally.
make.branches.geosse <- function(cache, control)
  make.branches.dtlik(cache$info, control)

derivs.geosse <- function(t, y, pars) {
  E_1 <- y[1]
  E_2 <- y[2]
  E_3 <- y[3]
  D_N1 <- y[4]
  D_N2 <- y[5]
  D_N3 <- y[6]

  sA  <- pars[1]
  sB  <- pars[2]
  sAB <- pars[3]
  xA  <- pars[4]
  xB  <- pars[5]
  dA  <- pars[6]
  dB  <- pars[7]

  ydot <- numeric(6)

  ## dE_1 / dt
  ydot[1] <- (-(sA + sB + xA + xB + sAB) * E_1 
              + xA * E_3 + xB * E_2 
              + sA * E_1 * E_2 + sB * E_1 * E_3 + sAB * E_2 * E_3)

  ## dE_2 / dt
  ydot[2] <- (-(sA + dA + xA) * E_2 
              + xA + dA * E_1 + sA * E_2 * E_2)

  ## E_3 / dt
  ydot[3] <- (-(sB + dB + xB) * E_3 
              + xB + dB * E_1 + sB * E_3 * E_3)

  ## dD_N1 / dt
  ydot[4] <- (-(sA + sB + sAB + xA + xB) * D_N1 
              + xA * D_N3 + xB * D_N2 
              + sA * (E_2 * D_N1 + E_1 * D_N2) 
              + sB * (E_3 * D_N1 + E_1 * D_N3)
              + sAB * (E_2 * D_N3 + E_3 * D_N2))

  ## dD_N2 / dt
  ydot[5] <- (-(sA + dA + xA) * D_N2 
              + dA * D_N1 + 2 * sA * D_N2 * E_2)

  ## dD_N3 / dt
  ydot[6] <- (-(sB + dB + xB) * D_N3 
              + dB * D_N1 + 2 * sB * D_N3 * E_3)

  ydot
}
