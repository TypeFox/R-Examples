make.all.branches.ou.pruning <- function(cache, control) {
  if (!is.ultrametric(cache$info$phy))
    warning("I am not sure that OU via pruning on non-ultrametric trees is calculated correctly")
  if (control$backend == "R") {
    branches.ou <- make.branches.ou(cache, control)
    function(pars, intermediates, preset=NULL)
      all.branches.matrix(pars, cache,
                          initial.conditions.bm.pruning,
                          branches.ou, preset)
  } else {
    # Force explicit selection of with- or without optimum functions.
    cache$info$name <- if (cache$with.optimum) "ou_opt" else "ou_noopt"
    make.all.branches.continuous(cache, control)
  }
}

## In ou, for a branch that ends with mean m (y[1]), the mean moves to
##   exp(len * alpha) * (m - theta) + theta
## And the variance becomes
##   (exp(2 * len * alpha) - 1) / (2 * alpha) + exp(2 * len * alpha) * v
## With normalising constant
##   exp(t * alpha)

## The last line of the second option comes from the limit of
##   (exp(2*len*alpha) - 1) * sigma2 / (2*alpha)
## as alpha -> 0 being len * sigma2
##   if ( alpha > 0 )
##     c(len * alpha + z,
##       exp(len * alpha) * (m - theta) + theta,
##       (exp(2*len*alpha) - 1) * sigma2 / (2*alpha) + exp(2*len*alpha) * v)
##   else
##     c(len * alpha + z,
##       exp(len * alpha) * (m - theta) + theta,
##       len * sigma2 + exp(2*len*alpha) * v)

## Here is the series of steps to get the "correct" equation for the
## no-optimim case, starting from rescale.R's equation
## t1 <- t0 + len # Branch base, backwards in time
## t.end <- t.max - t0
## t.start <- t.end - len
## cmp1 <-
##   exp(-2 * alpha * t.max) * (
##     expm1(2 * alpha * t.end) - expm1(2 * alpha * t.start)) /
##       (2 * alpha)
## cmp2 <- exp(-2 * alpha * (t.max - t.start)) *
##   expm1(2 * alpha * len)/(2*alpha)
## cmp3 <- exp(-2 * alpha * (t0 + len)) *
##   expm1(2 * alpha * len)/(2*alpha)
## cmp4 <- -exp(-2 * alpha * t0) *
##   expm1(-2 * alpha * len)/(2*alpha)
make.branches.ou <- function(cache, control) {
  with.optimum <- cache$with.optimum

  function(y, len, pars, t0, idx) {
    m <- y[[1]]
    v <- y[[2]]
    z <- y[[3]]

    sigma2 <- pars[[1]]
    alpha  <- pars[[2]]

    if (alpha == 0) { # see branches.bm.pruning()
      list(z, c(m, v + len * sigma2, 0))
    } else if (with.optimum) {
      theta  <- pars[[3]]
      list(len * alpha + z,
           c(exp(len * alpha) * (m - theta) + theta,
             expm1(2*len*alpha) * sigma2 / (2*alpha) +
             exp(2*len*alpha) * v,
             0))
    } else { # no optimum, based on scaling
      len.scaled <- exp(-2 * alpha * t0) * expm1(-2 * alpha * len) /
        (2*alpha)
      list(z,
           c(m, v - sigma2 * len.scaled, 0))
    }
  }
}
