## ------------------------------------------------------------------------
# Part 1:  The likelihood function
triangular.like <- function(b, dist, w.lo, w.hi, series="", expansions=0, scale=TRUE){
  L <- (2/b)*(1 - dist/b)
  L[ L < 0 ] <- 0
  L
  }


## ------------------------------------------------------------------------
# Part 2:  The starting values, limits, and names of parameters in the likelihood
triangular.start.limits <- function(dist, expansions, w.lo, w.hi){
  list(start=max(dist)*.75,
       lowlimit=w.lo,
       highlimit=w.hi,
       names="Max")
  }


## ---- fig.width=6, fig.height=4------------------------------------------
# A function to generate triangular random deviates
rtriang <- function(n, b){
  x <- seq(0, b, length=500)
  CDF <- 2*x/b - (x/b)^2
  u <- runif(n)
  r <- approx( CDF, x, xout=u )$y
  }

# Simulated vector of distances
set.seed(123)
d <- rtriang(500, 100)  # true b = 100
hist(d)

## ---- fig.width=6, fig.height=4------------------------------------------
# Fit detection function with user-defined "triangular" likelihood
# Requires the F.dfunc.estim function from Rdistance
require(Rdistance)
tri.dfunc <- F.dfunc.estim( d, likelihood="triangular", w.hi=150 )
plot(tri.dfunc)

## ------------------------------------------------------------------------
tri.dfunc$g.x.scl*tri.dfunc$param / 2

## ------------------------------------------------------------------------

ESW(tri.dfunc)

