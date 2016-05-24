### This file is for Cov(z_nk) and Cov(z_nk - z_n1).

### Partial logL. Return a matrix with dimension M * N.
# partial.q <- function(x, PI, MU, S, post.z, logit.PI = TRUE)

### Compute posterior. Return a matrix with dimension N * K.
# postPI <- function(x, emobj)

### Partial all z_nk. Return a list of length nrow(x).
partial.post.z <- function(x, PI, MU, S, post.z, index.N){
  invS <- apply(S, 3, solve)
  dim(invS) <- dim(S)

  ret <- list()
  for(i.n in 1:length(index.N)){
    ret[[i.n]] <- partial.post.z.1(x, PI, MU, S, post.z, invS, index.N[i.n])
  }

  ret
} # End of partial.post.z().

### Partial one z_nk. Return a matrix with dimension K * M for index i.n.
partial.post.z.1 <- function(x, PI, MU, S, post.z, invS, i.n){
  K <- length(PI)

  ### Make a K * K weight matrix.
  Zbuf.w <- list()
  for(i.k in 1:K){
    tmp <- rep(-post.z[i.n, i.k], K)
    tmp[i.k] <- tmp[i.k] + 1
    Zbuf.w[[i.k]] <- post.z[i.n,] / PI * tmp
  }

  ### For PI
  Zbuf.PI <- list()
  if(K > 1){
    for(i.k in 1:(K - 1)){
      Zbuf.PI[[i.k]] <- Zbuf.w[[i.k]][-K] +
                         post.z[i.n, K] / PI[K] * post.z[i.n, i.k]
    }
    Zbuf.PI[[K]] <- -colSums(do.call("rbind", Zbuf.PI))
  } else{
    Zbuf.PI[[1]] <- vector(mode = "numeric", length = 0)
  }

  ### For MU and S.
  Zbuf.MU <- list()
  Zbuf.S <- list()
  for(i.k in 1:K){
    tmp.mu <- list()
    tmp.s <- list()
    for(j.k in 1:K){
      tmp <- partial.f.mu.s(x, MU, S, invS, i.n, j.k)
      tmp.mu[[j.k]] <- Zbuf.w[[i.k]][j.k] * tmp$mu
      tmp.s[[j.k]] <- Zbuf.w[[i.k]][j.k] * tmp$s
    }
    Zbuf.MU[[i.k]] <- do.call("c", tmp.mu)
    Zbuf.S[[i.k]] <- do.call("c", tmp.s)
  }

  ### For Zbuf. Combine all PI, MU, and S.
  Zbuf <- list()
  for(i.k in 1:K){
    Zbuf[[i.k]] <- c(Zbuf.PI[[i.k]], Zbuf.MU[[i.k]], Zbuf.S[[i.k]])
  }
  Zbuf <- do.call("rbind", Zbuf)

  Zbuf
} # End of partial.post.z.1().

partial.f.mu.s <- function(x, MU, S, invS, i.n, j.k){
  f <- dmvn(x[i.n,], MU[j.k,], var2LTsigma(S[,, j.k]))
  tmp.x.minus.mu <- matrix(x[i.n,] - MU[j.k,], ncol = 1)
  tmp.mu <- invS[,, j.k] %*% tmp.x.minus.mu
  # tmp.s <- invS[,, j.k] %*% (tmp.x.minus.mu %*% t(tmp.x.minus.mu) %*%
  #                            invS[,, j.k] - diag(1, nrow(MU)))
  # tmp.s <- var2LTsigma(tmp.s) / 2
  tmp.s <- var2LTsigma(tmp.mu %*% t(tmp.mu) - invS[,, j.k]) / 2
  list(mu = f * tmp.mu, s = f * tmp.s)
} # End of partial.f.mu.s().

get.cov.param <- function(x, emobj, post.z){
  K <- emobj$nclass
  PI <- emobj$pi
  MU <- emobj$Mu
  S <- LTSigma2variance(emobj$LTSigma)

  ### cov of param = {PI, MU, S}.
  pl <- partial.q(x, PI, MU, S, post.z, logit.PI = FALSE)

  ### Get cov of param
  I <- pl %*% t(pl)
  cov <- ginv(I)

  list(I = I, cov = cov)
} # End of get.cov.param().

get.cov.post.z <- function(x, emobj, post.z, cov.param = NULL){
  K <- emobj$nclass
  PI <- emobj$pi
  MU <- emobj$Mu
  S <- LTSigma2variance(emobj$LTSigma)

  ### nabla of post.z.
  index.N <- 1:nrow(x)
  nabla.post.z <- partial.post.z(x, PI, MU, S, post.z, index.N)

  ### Get cov of post.z
  if(is.null(cov.param)){
    cov.param <- get.cov.param(x, emobj)$cov
  }
  cov.post.z <- lapply(nabla.post.z, function(x) x %*% cov.param %*% t(x))
  cov.post.z
} # End of get.cov.post.z().

get.cov.logit.z <- function(x, emobj, post.z, cov.param = NULL,
    cov.post.z = NULL){
  K <- emobj$nclass

  ### nabla of logit.z.
  index.N <- 1:nrow(x)
  nabla.logit.z <- lapply(index.N, function(i){ partial.logit.p(post.z[i,]) })

  ### Get cov of logit.z
  if(is.null(cov.param)){
    cov.param <- get.cov.param(x, emobj, post.z)$cov
  }
  if(is.null(cov.post.z)){
    cov.post.z <- get.cov.post.z(x, emobj, post.z, cov.param = cov.param)
  }
  if(length(nabla.logit.z) != length(cov.post.z)){
    stop("nabla.logit.z and cov.post.z are not of the same length.")
  }
  cov.logit.z <- lapply(1:length(nabla.logit.z),
                        function(i){
                          nabla.logit.z[[i]] %*%
                          cov.post.z[[i]][1:(K-1), 1:(K-1)] %*%
                          t(nabla.logit.z[[i]])
                        })
  cov.logit.z
} # End of get.cov.logit.z().

get.logor.stat <- function(x, emobj, post.z, cov.param = NULL,
    cov.post.z = NULL, cov.logit.z = NULL){
  K <- emobj$nclass

  ### Get cov of logor
  if(is.null(cov.param)){
    cov.param <- get.cov.param(x, emobj, post.z)$cov
  }
  if(is.null(cov.post.z)){
    cov.post.z <- get.cov.post.z(x, emobj, post.z, cov.param = cov.param)
  }
  if(is.null(cov.logit.z)){
    cov.logit.z <- get.cov.logit.z(x, emobj, post.z, cov.param = cov.param,
                                   cov.post.z = cov.post.z)
  }

  ### Get log odds ratio and its cov matrix.
  index.N <- 1:nrow(x)
  logit.p <- log(post.z[index.N,] / (1 - post.z[index.N,]))
  A <- cbind(rep(-1, K - 1), diag(1, K - 1))
  log.or <- logit.p %*% t(A)
  cov.log.or <- lapply(cov.logit.z, function(x) A %*% x %*% t(A))

  ### Get statistics and p-values.
  df <- rep(NA, nrow(log.or))
  test.stat <- rep(NA, nrow(log.or))
  for(i in 1:nrow(log.or)){
    ts <- matrix(log.or[i,] - 1, nrow = 1)    # check if log.or == 1
    if(all(is.finite(ts)) && all(is.finite(cov.log.or[[i]]))){
      # df[i] <- sum(eigen(cov.log.or[[i]], only.values = TRUE)$values > 0)
      df[i] <- as.integer(rankMatrix(cov.log.or[[i]]))
      test.stat[i] <- ts %*% ginv(cov.log.or[[i]]) %*% t(ts)
    }
  }

  ### Return.
  list(log.or = log.or, cov.log.or = cov.log.or, df = df,
       test.stat = test.stat)
} # End of get.logor.stat().

