# From SamplerCompare, (c) 2010 Madeleine Thompson

# Covariance-matching slice sampler, see online help and the following
# technical report for more information:
#
#   Thompson, M. B. and Neal, R. M. (2010). Covariance-adaptive
#   slice sampling. Technical Report TR-1002, Dept. of Statistics,
#   University of Toronto.

cov.match.step <- function(target.dist, x0, y0, tuning=1, theta=0.95) {
  p <- length(x0)
  evals <- 0
  grads <- 0
  y.max <- -Inf                       # max. log density so far this obs.
  if (y0>y.max)
    y.max <- y0
  y <- y0 - rexp(1)                   # slice level
  crumb.R <- diag(1/tuning, nrow=p)   # Cholesky factor of crumb precision
  nprop <- 0                          # (F_k in Thompson & Neal)

  # Repeatedly draw crumbs and proposals till a proposal is accepted.

  repeat {
    nprop <- nprop + 1

    # Draw a crumb.

    crumb.std <- rnorm(p)                           # z_1 in Thompson & Neal
    crumb.offset <- backsolve(crumb.R, crumb.std)
    crumb <- crumb.offset + x0                      # c_k in Thompson & Neal

    # Update posterior.R, the Cholesky factor of the precision
    # of the posterior distribution for x0, given that we've observed
    # nprop crumbs.  This corresponds to eqn. 28 in Thompson &
    # Neal.  Then, update the standardized posterior mean (\bar c_k^*)
    # and posterior mean (\bar c_k) using eqns. 31 and 32 of
    # Thompson & Neal.

    if (nprop==1) {
      posterior.std <- crumb.std
      posterior.mean <- crumb.offset
      posterior.R <- crumb.R
    } else {
      posterior.std <- (t(crumb.R) %*% ( crumb.R %*% crumb.offset ) +
                        t(posterior.R) %*% ( posterior.R %*% posterior.mean))
      posterior.R <- chud(sqrt(1+theta)*posterior.R, sqrt(alpha)*g)
      posterior.mean <- backsolve(posterior.R,
        backsolve(posterior.R, posterior.std, transpose=TRUE))
    }

    # Draw a proposal, x1, from the posterior for x0.
    # If it is inside the slice, accept it.

    Delta <- backsolve(posterior.R, rnorm(p)) 
    x1 <- Delta + posterior.mean + x0          # proposal
    evals <- evals + 1
    y1 <- target.dist$log.density(x1)          # log-density at proposal
    if (y1>y.max)
      y.max <- y1
    if (y1>=y)                                 # inside the slice?
      break

    # The proposal is not in the slice.  Compute the gradient at
    # the proposal and try to adapt the crumb covariance so that
    # it will produce a more suitable posterior distribution.

    grads <- grads + 1
    g <- array(target.dist$grad.log.density(x1), c(p,1))  # gradient
    if (all(is.finite(g))) {
      x2 <- x1 + g/twonorm(g) * twonorm(x1-crumb)         # arbitrary point
      evals <- evals + 1
      y2 <- target.dist$log.density(x2)                   # log density at x2
      if (!is.na(y2) && y2>y.max)
        y.max <- y2
      d2 <- twonorm(x1-crumb)
      curv <- -2/d2^2 * ( y2-y1-twonorm(g)*d2 )           # kappa, eqn. 14
      cut.max <- -1/2*twonorm(g)^2/curv + twonorm(g)^2/curv + y1
      if (!is.na(cut.max) && cut.max > y.max)
        y.max <- cut.max
      gg <- sum(g*g)
      Rg <- posterior.R %*% g
      alpha <- (1.5/(y.max-y) * curv - (1+theta)*sum(Rg*Rg)/gg)/gg  # eqn. 23
    } else {
      alpha <- 0
    }

    # If alpha could not be computed due to numeric reasons, zero
    # it (and the gradient) out before updating the crumb precision.

    if (!isTRUE(alpha>0) || !all(is.finite(sqrt(alpha)*g))) {
      g <- array(0, c(p,1))
      alpha <- 0
    }

    # Update the crumb precision.  Perform numeric stability
    # checks before committing to it.

    crumb.R.new <- chud(sqrt(theta)*posterior.R, sqrt(alpha)*g)

    if (!all(is.finite(crumb.R.new)) ||
        is.nan(k <- rcond(crumb.R.new)) ||
        k<0.00001) {
      alpha <- 0
    } else {
      crumb.R <- crumb.R.new
    }
  }
  return(list(x=x1, y=y1, evals=evals, grads=grads))
}

# compounded.sampler is not used directly because the signature of
# its return value has ..., causing spurious warnings if the .Rd file
# for cov.match.sample has the complete signature.

cov.match.sample <- function(target.dist, x0, sample.size, tuning=1, theta=1,
                             limit=length(x0)*100) {
  s <- compounded.sampler(list(cov.match.step), 'Cov. Matching (inner)')
  s(target.dist, x0, sample.size, tuning=tuning, theta=theta, limit=limit)
}
attr(cov.match.sample, 'name') <- 'Cov. Matching'
