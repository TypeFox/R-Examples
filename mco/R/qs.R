##
## qs.R - Quality measures 
##
## Authors:
##  Heike Trautmann  <trautmann@statistik.uni-dortmund.de>
##  Detlef Steuer    <detlef.steuer@hsu-hamburg.de>
##  Olaf Mersmann    <olafm@statistik.uni-dortmund.de>
##

## Helper functions first:

## Normalize front so that all dimension are in [0, 1]
normalizeFront <- function(front, minval, maxval) {
  if (missing(minval))
    minval <- apply(front, 2, min)
  if (missing(maxval))
    maxval <- apply(front, 2, max)
  ## FIXME: This is ugly!
  t((t(front) - minval)/(maxval - minval))
}

## Squared distance
distance2 <- function(x, y)
  crossprod((x - y))

## Squared distance to front
distanceToFront2 <- function(x, front)
  min(sapply(1:nrow(front), function(i) { crossprod(x - front[i,]) }))

generationalDistance <- function(x, o) {
  front <- paretoFront(x)
  truefront <- paretoFront(o)

  ## Normalize front:
  maxval <- apply(truefront, 2, max)
  minval <- apply(truefront, 2, min)

  nfront <- normalizeFront(front, minval, maxval)
  ntruefront <- normalizeFront(truefront, minval, maxval)

  ## Calculate criterion:
  d <- sapply(1:nrow(nfront), function(i) distanceToFront2(nfront[i,], ntruefront))
  return(sqrt(sum(d))/nrow(nfront))
}

generalizedSpread <- function(x, o) {
  front <- paretoFront(x)
  truefront <- paretoFront(o)
  
  ## Normalize front:
  maxval <- apply(truefront, 2, max)
  minval <- apply(truefront, 2, min)

  nfront <- normalizeFront(front, minval, maxval)
  ntruefront <- normalizeFront(truefront, minval, maxval)

  K <- nrow(nfront)
  N <- nrow(ntruefront)
  ## Calculate extreme values:
  nobj <- ncol(front)
  
##  extreme <- matrix(0, ncol=nobj, nrow=nobj)
##  for (i in 1:nobj) {
##    o <- order(ntruefront[,i])
##    for (j in 1:nobj) {
##      extreme[i,j] <- ntruefront[o,][N, j]
##    }
##  }
  extreme <- sapply(1:nobj, function(i) ntruefront[which.max(ntruefront[,i]),])

  ## Lexographically sort front:
  for (i in nobj:1)
    nfront <- nfront[order(nfront[,i]),]
  
  if (distance2(nfront[1,], nfront[K,]) == 0) {
    return (0.0)
  } else {
    dmean <- mean(sapply(1:K, function(i) sqrt(distanceToFront2(nfront[i,], nfront[-i,]))))
    dextr <- sum(sapply(1:nobj, function(i) sqrt(distanceToFront2(extreme[i,], nfront))))
    mean <-  sum(sapply(1:K, function(i) sqrt(distanceToFront2(nfront[i,], nfront[-i,]))-dmean))
    return ((dextr + mean)/(dextr + K*dmean))                 
  }
}

dominatedHypervolume <- function(x, ref) {
  ## Extract pareto front:
  front <- paretoFront(x)
  
  ## Possibly infer reference point:
  if (missing(ref))
    ref <- apply(front, 2, max)

  ## Sanity checks:
  if (!is.matrix(front))
    stop("Pareto front must be a matrix")
  if (ncol(front) != length(ref))
    stop("Reference point and front must have the same dimension.")

  ## Note the transopse. do_hv() needs the front in row major format.
  .Call(do_hv, t(front), ref)
}

epsilonIndicator <- function(x, o) {
  x.front <- paretoFront(x)
  o.front <- paretoFront(o)

  if (!is.matrix(x.front))
    stop("'x' must be a valid pareto front.")
  if (!is.matrix(o.front))
    stop("'o' must be a valid pareto front.")#

  if (any(x.front < 0) || any(o.front < 0))
    stop("Epsilon Indicator only works for fronts which are strictly positive.")
  
  .Call(do_eps_ind, x.front, o.front)
}
