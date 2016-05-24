## functions for approximate Bayes factors

## wakefield ABF inverted to be BF for alternate model relative to null
## assumes prior mean is zero
abf.Wakefield <- function(beta, se, priorsd) return(sqrt(se^2/(se^2 + priorsd^2)) * exp((beta/se)^2/2 * priorsd^2/(se^2+priorsd^2)))

## code used for 2011 AJHG paper used fixed quadrature grid, xx <- seq(from = -0.5, to = 0.5, length.out = 1001)

## ABF for normal prior calculated numerically, should be equal to that calculated by abf.wakefield()
abf.normal <- function(beta, se, priorscale, gridrange = 3, griddensity = 20) {
  gridrange <- max(1, gridrange)
  griddensity <- max(1, griddensity)
  xlim <- range(c(beta - gridrange*se, beta + gridrange*se, -gridrange*priorscale, gridrange*priorscale), na.rm = TRUE)
  xx <- seq(from = xlim[1], to = xlim[2], length.out = max(3, ceiling(griddensity*(xlim[2] - xlim[1])/min(se, na.rm = TRUE))))
  dxx <- xx[2] - xx[1]
  return(sapply(1:length(beta), function(idx) {
    sum(dxx * dnorm(xx/priorscale)/priorscale * dnorm(beta[idx], xx, se[idx]))/dnorm(beta[idx], 0, se[idx])
  }))
}

## ABF for t distribution prior (a priori beta/priorscale ~ standard t distribution)
abf.t <- function(beta, se, priorscale, df = 1, gridrange = 3, griddensity = 20) {
  gridrange <- max(1, gridrange)
  griddensity <- max(1, griddensity)
  xlim <- range(c(beta - gridrange*se, beta + gridrange*se, -gridrange*priorscale, gridrange*priorscale), na.rm = TRUE)
  xx <- seq(from = xlim[1], to = xlim[2], length.out = max(3, ceiling(griddensity*(xlim[2] - xlim[1])/min(se, na.rm = TRUE))))
  dxx <- xx[2] - xx[1]
  return(sapply(1:length(beta), function(idx) {
    sum(dxx * dt(xx/priorscale, df = df)/priorscale * dnorm(beta[idx], xx, se[idx]))/dnorm(beta[idx], 0, se[idx])
  }))
}
