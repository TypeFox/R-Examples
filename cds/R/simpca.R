#' Simulate Data with a Specific Principal Components Structure and Response
#' Style Contamination
#' 
#' Simulate normally distributed data with specific covariance structure and
#' randomly sampled means. Adds response style contamination.
#' @param nr.indv Numeric vector of group sizes.
#' @param m Integer; then number of variables to simulate.
#' @param q Integer; the rating scale used \code{1:q}.
#' @param R List with entry named 'R' which is the simulated correlation matrix
#' @param err.coeff Standard error for each variable, added unto \code{R}.
#' @param alphamat Matrix containing splines coefficients for te construction
#' of respone styles.
#' @param randomize logical; should the rows of the data be randomly permuted
#' or not?
#' @export simpca
simpca <- function(nr.indv = rep(200, 5), m = 10, q = 7, R = rcormat(m = m), err.coeff = 0.1,
                    alphamat = rbind(c(0.5, 2, 4), c(10, 2, 10), c(1, 2, 1), 
                    c(4, 2, 0.5), c(0.1, 2, 0.1))[1:length(nr.indv),], randomize = FALSE){
  Rmat <- R$R
  Sigma <- Rmat*matrix(err.coeff^2, nrow = m, ncol = m)
  mns <- runif(m)
#   struct <- eigen(Sigma)
  n <- sum(nr.indv)
  if(nrow(alphamat) != length(nr.indv)) stop('Inconsistent arguments: nr.indv and alphamat imply different K')
  K <- length(nr.indv)
  Xstar <- MASS::mvrnorm(n = n, mu = mns, Sigma = Sigma)
#   Xstarnull <- sweep(Xstar, MARGIN = 2, STATS = mns, FUN = '-')
#   Xstd <- sweep(Xstarnull, MARGIN = 2, STATS = sqrt(diag(Sigma)), FUN = '/')
  grp.ind <- rep(1:K, times = nr.indv)
  tvec <- quantile(Xstar, probs = c(0.01, 0.99))
  dif <- diff(range(tvec))
  tvec <- c(tvec[1], mean(tvec), tvec[2])
  rangeX <- range(Xstar)
  rangeHi <- max(abs(rangeX))
  cpoints <- apply(alphamat, 1, create.rs, nr.scale = q, tvec = tvec, 
                    xvec = seq(from = tvec[1], to = tvec[3], length = q + 1))
  cpoints <- cpoints*dif + tvec[1]
  cpoints[1,] <- -Inf
  cpoints[q + 1,] <- Inf
  Xsplit <- split(as.data.frame(Xstar), f = grp.ind) 
  for(i in 1:K){
      Xsplit[[i]] <- matrix(as.numeric(cut(unlist(Xsplit[[i]]), breaks = cpoints[,i], labels = 1:q)), 
                            nrow = nr.indv[i], ncol = m)
  }
  X <- do.call(rbind, Xsplit)  
  if(randomize){
    ord <- sample(1:n, n)
    X <- X[ord, ]
    grp.ind <- grp.ind[ord]
    Xstar <- X[ord, ]
  }
  out <- list(sim = addbounds(X, q = q), mu = mns, resp.style.alpha = alphamat, 
              grp.ind = grp.ind, G = indmat(grp.ind, K = K), m = m, 
              scales = 1:q, data = X, pre.rs = Xstar, true.tau = NULL, Sigma = Sigma, err.coeff = err.coeff, R = R)
  class(out) <- c("dspca", "dsdata")
  
#   pts <- seq(from = -rangeHi, to = rangeHi, length = 500)
#   plot(0, 0, type = "n", xlim = range(Xstar), ylim = c(0, dnorm(0, sd = err.coeff)))
#   for(i in 1L:m)
#     lines(pts, dnorm(pts, mean=mns[i], sd=sqrt(Sigma[i,i])))
#   for(k in 1:K)
#     abline(v = cpoints[,k], col = k)
  out
}
