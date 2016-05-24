### This file contains files for estimating parameters emperically.


### Estimate MU based on classified data.
estimate.MU <- function(X.spmd, CLASS.spmd, K){
  p <- ncol(X.spmd)
  MU.CLASS <- matrix(0, nrow = p, ncol = K)

  for(i.k in 1:K){
    tmp.id <- which(CLASS.spmd == i.k)
    tmp.n.id <- spmd.allreduce.double(as.double(length(tmp.id)),
                                      double(1), op = "sum")

    if(length(tmp.id) > 0){
      tmp.MU <- colSums(X.spmd[tmp.id,])
    } else{
      tmp.MU <- rep(0.0, p)
    }
    tmp.MU <- spmd.allreduce.double(tmp.MU, double(p), op = "sum")

    MU.CLASS[, i.k] <- tmp.MU / tmp.n.id
  }

  MU.CLASS
} # End of estimate.MU().


### Estimate SIGMA based on classified data.
#my.estimate.sigma <- function(n, k, X.spmd, MU){
#  x <- matrix(X.spmd[n,] - MU[, k], nrow = 1)
#  as.vector((t(x) %*% x))
#} # End of my.estimate.sigma().

estimate.SIGMA <- function(X.spmd, MU, CLASS.spmd, K){
  p <- ncol(X.spmd)
  p.2 <- p * p
  SIGMA.CLASS <- NULL

  for(i.k in 1:K){  
    tmp.id <- which(CLASS.spmd == i.k)
    tmp.n.id <- spmd.allreduce.double(as.double(length(tmp.id)),
                                      double(1), op = "sum")

    if(length(tmp.id) == 1){
      tmp.X.spmd <- X.spmd[tmp.id,] - MU[, i.k]
      dim(tmp.X.spmd) <- c(1, p)
      tmp.SIGMA <- crossprod(tmp.X.spmd)
    } else if(length(tmp.id) > 1){
# Version 1:
#      tmp.SIGMA <- rowSums(do.call("cbind",
#                                   lapply(tmp.id, my.estimate.sigma,
#                                          i.k, X.spmd, MU)))
# Version 2:
#      tmp.X.spmd <- t.X.spmd[tmp.id,] - MU[, i.k]
#      tmp.SIGMA <- tmp.X.spmd %*% t(tmp.X.spmd)
# Version 3:
      tmp.X.spmd <- W.plus.y(X.spmd[tmp.id,], -MU[, i.k],
                               length(tmp.id), p)
      tmp.SIGMA <- crossprod(tmp.X.spmd)
    } else{
      tmp.SIGMA <- rep(0.0, p.2)
    }
    tmp.SIGMA <- spmd.allreduce.double(tmp.SIGMA, double(p.2), op = "sum")

    SIGMA.CLASS[[i.k]] <- matrix(tmp.SIGMA / (tmp.n.id - 1), ncol = p)
  }

  SIGMA.CLASS
} # End of estimate.SIGMA().


### This function collects N.CLASS
get.N.CLASS <- function(K){
  tmp.n.class <- tabulate(.pmclustEnv$CLASS.spmd, nbins = K)
  spmd.allreduce.integer(as.integer(tmp.n.class), integer(K), op = "sum")
} # End of get.N.CLASS().


### Copy from dmvnorm() in mvtnorm package.
#logdmvnorm <- function(x, mean, sigma){
#  x <- matrix(x, nrow = 1)
#  distval <- mahalanobis(x, center = mean, cov = sigma)
#  logdet <- sum(log(eigen(sigma, symmetric = TRUE, only.values = TRUE)$values))
#  logretval <- -(ncol(x) * log(2 * pi) + logdet + distval)/2
#  logretval
#} # End of logdmvnorm()

