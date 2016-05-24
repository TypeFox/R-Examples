# this function was copied from the package mvtnorm to reduce the
# number of packages that birdring is dependent on.
#-------------------------------------------------------
# please cite the following reference: 
#Alan Genz, Frank Bretz, Tetsuhisa Miwa, Xuefei Mi,
#Friedrich Leisch, Fabian Scheipl, Torsten Hothorn
#(2013). mvtnorm: Multivariate Normal and t
#Distributions. R package version 0.9-9995. URL
#http://CRAN.R-project.org/package=mvtnorm

#Alan Genz, Frank Bretz (2009), Computation of
#Multivariate Normal and t Probabilities. Lecture Notes
#in Statistics, Vol. 195., Springer-Verlage, Heidelberg.
#ISBN 978-3-642-01688-2
#------------------------------------------------

dmvnorm <- function (x, mean, sigma, log = FALSE) 
{
  if (is.vector(x)) {
    x <- matrix(x, ncol = length(x))
  }
  if (missing(mean)) {
    mean <- rep(0, length = ncol(x))
  }
  if (missing(sigma)) {
    sigma <- diag(ncol(x))
  }
  if (NCOL(x) != NCOL(sigma)) {
    stop("x and sigma have non-conforming size")
  }
  if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps), 
                   check.attributes = FALSE)) {
    stop("sigma must be a symmetric matrix")
  }
  if (length(mean) != NROW(sigma)) {
    stop("mean and sigma have non-conforming size")
  }
  distval <- mahalanobis(x, center = mean, cov = sigma)
  logdet <- sum(log(eigen(sigma, symmetric = TRUE, only.values = TRUE)$values))
  logretval <- -(ncol(x) * log(2 * pi) + logdet + distval)/2
  if (log) 
    return(logretval)
  exp(logretval)
}