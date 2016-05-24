#'
#' Transform back multiple regression coefficients to unscaled regression coefficients
#' Original question posed by Mark Seeto on the R mailing list.
#' 
#'@author M.Suzen
#'@note 2015-04-10
#'@param X, unscaled design matrix without the intercept, m by n matrix
#'@param Y, unscaled response, m by 1 matrix
#'@param betas.scaled, coefficients vector of multiple regression, first term is the intercept
#'@examples
#'  set.seed(4242)
#'  X            <- matrix(rnorm(12), 4, 3)
#'  Y            <- matrix(rnorm(4), 4, 1)
#'  betas.scaled <- matrix(rnorm(3), 3, 1)
#'  betas        <- scaleBack.lm(X, Y, betas.scaled)
scaleBack.lm <- function(X, Y, betas.scaled) {
  numB     <- length(betas.scaled)
  betas    <- rep(0.0, numB)
  X.s      <- attributes(scale(X))
  muX      <- X.s$`scaled:center`
  sigmaX   <- X.s$`scaled:scale`
  Y.s      <- attributes(scale(Y))
  muY      <- Y.s$`scaled:center`
  sigmaY   <- Y.s$`scaled:scale`
  vect     <- sapply(1:(numB-1), function(inx) muX[inx]*sigmaY * betas.scaled[inx+1]/sigmaX[inx])
  betas             <- matrix(0.0, numB, 1)
  betas[1]          <- betas.scaled[1]*sigmaY+muY-sum(vect) # intercept
  betas[2:numB]     <- sapply(1:(numB-1), function(inx) betas.scaled[inx+1]*sigmaY/sigmaX[inx])
  return(betas)
}
