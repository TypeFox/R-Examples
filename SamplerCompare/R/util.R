# From SamplerCompare, (c) 2010 Madeleine Thompson

# util.R is a place where assorted utility functions find a home.

# Not currently used; draws n rows uniformly from the unit sphere in R^p.
# Read http://en.wikipedia.org/wiki/Hypersphere for more.

rsphere <- function(n,p) {
  r <- runif(n)^(1/p)
  x <- array(rnorm(n*p), c(n,p))
  x <- x / sqrt(rowSums(x^2)) * r
  return(x)
}

# Sooner or later, every program I write gets a logsumexp function.
# Sometimes two.

logsumexp <- function(x) {
  m <- max(x)
  if (is.na(m) || m==-Inf)
    return(m)
  return(m+log(sum(exp(x-m))))
}

# Returns the two-norm of the argument vector.  Another function I
# use all the time.

twonorm <- function(x) sqrt(sum(x^2))
