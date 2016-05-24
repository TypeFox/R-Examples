## Helper function to handle (1 + x*shape)^(-1/shape) as shape -> 0
nzsh <- function(x, shape) {
  ifelse(shape == 0, exp(-x), exp((-1/shape) * log1p(x * shape)))
}


## Function to help deal with design matrix
adjScale <- function(x) {
  truemeans <- as.numeric(colMeans(x))
  truevars <- as.numeric(apply(x, 2, sd))
  adjmeans <- ifelse(truevars == 0, 0, truemeans)
  adjvars <- ifelse(truevars == 0, truemeans, truevars)
  if(ncol(x) == 1)
    adjmeans <- 0
  x <- t((t(x) - adjmeans) / adjvars)
  out <- list(x, truemeans, truevars, adjmeans, adjvars)
  names(out) <- c("mat", "truemeans", "truevars", "adjmeans", "adjvars")
  out
}


## Returns expected inverse fisher information matrix for GEVr or Gumbel distribution.
gevrFisher <- function(data, theta, gumbel = FALSE) {
  data <- as.matrix(data)
  R <- ncol(data)
  N <- nrow(data)
  loc <- theta[1]
  scale <- theta[2]
  if(!gumbel) {
    shape <- theta[3]
    gr1 <- gamma(R + shape + 1) / gamma(R)
    gr2 <- gamma(R + 2*shape + 1) / gamma(R)
    A <- (((1+shape)^2)*gr2)/((scale^2)*(1+2*shape))
    B <- (-1/((scale^2)*shape*(1+2*shape)))*(((1+shape)^2)*gr2 - (1+2*shape)*gr1)
    C <- (1/(scale*(shape^2)*(1+2*shape)))*(((1+2*shape)*gr1)*(shape*digamma(R+shape+1) + (shape^2+shape+1)/(1+shape)) - ((1+shape)^2)*gr2)
    D <- (1/((scale^2)*(shape^2)*(1+2*shape)))*(R*(1+2*shape) - 2*(1+2*shape)*gr1 + ((1+shape)^2)*gr2)
    E <- (1/(scale*((-shape)^3)*(1+2*shape)))*(R*(-shape)*(1+2*shape)*digamma(R+1) + (1+2*shape)*gr1*(shape*digamma(R+shape+1) + (1+(1+shape)^2)/(1+shape)) - ((1+shape)^2)*gr2 - R*(1+2*shape))
    F <- (1/(((-shape)^4)*(1+2*shape)))*((2*(1+2*shape)*gr1)*((-shape)*digamma(R+shape+1) - (shape^2+shape+1)/(1+shape)) + ((1+shape)^2)*gr2 + (R*(1+2*shape))*(1 + 2*shape*digamma(R+1) + (shape^2)*(1 + trigamma(R+1) + ((digamma(R+1))^2))))
    info <- solve(matrix(c(A, B, C, B, D, E, C, E, F), nrow=3, ncol=3)) / N
  } else {
    Br <- R * digamma(R + 1)
    Cr <- R * (digamma(R + 1)^2 + trigamma(R + 1) + 1)
    info <- (scale^2 / (N * (R * Cr - Br^2))) * matrix(c(Cr, Br, Br, R), nrow=2, ncol=2)
  }
  info
}


## Outputs matrix with row contributions to score. Need to sum over the columns to get full score.
gevrScore <- function(data, theta) {
  data <- as.matrix(data)
  R <- ncol(data)
  N <- nrow(data)
  loc <- theta[1]
  scale <- theta[2]
  shape <- theta[3]
  z <- (data - loc) / scale
  dLoc <- rowSums((((1/shape)+1)*shape) / (scale*(1+shape*z)))
  dLoc <- dLoc - nzsh(z[,R], shape) / (scale*(1+(shape*z[,R])))
  dScale <- rowSums((((1/shape)+1)*shape*z) / (scale*(1+shape*z)))
  dScale <- dScale - (R/scale) - z[,R] * nzsh(z[,R], shape) / (scale*(1+shape*z[,R]))
  dShape <-  rowSums(log(1+shape*z) / shape^2 - (((1/shape)+1)*z) / (1+shape*z))
  dShape <- dShape - nzsh(z[,R], shape) * ((log(1+shape*z[,R]) / shape^2) - z[,R] / (shape*(1+shape*z[,R])))
  score <- matrix(c(dLoc, dScale, dShape), ncol = 3)
  score
}


## Calculates the score test statistic
gevrTestStat <- function(z, information) {
  data <- z$data
  theta <- z$par.ests
  R <- ncol(data)
  N <- nrow(data)
  u <- gevrScore(data, theta)
  u <- colSums(u)
  if(information == "observed") {
    info <- z$varcov
  } else {
    info <- gevrFisher(data, theta)
  }
  stat <- as.vector((1/N) * t(u) %*% info %*% u)
  stat
}
