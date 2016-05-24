
# This function creates G matrix needed for taking
# derivatives of symmetric matrices
# Namely, vec(X) = G %*% vech(X)
Gmat <- function(p){
  M <- p * (p + 1) / 2
  G <- matrix(NA, ncol = M, nrow = p * p)
  n <- 1

  Z <- rep(0, M)
  for(a in 1:p){
    for(b in 1:p){
      Zn <- Z
      i1 <- max(a, b)      
      i2 <- min(a, b)
      ind <- M - (p - i2 + 1) * (p - i2 + 2) / 2 + i1 - i2 + 1
      Zn[ind] <- 1
      G[n,] <- Zn
      n <- n + 1
    }
  }
  return(G)
} # End of Gmat().


### partial.q. Return a matrix with dimension M * N.
partial.q <- function(x, PI, MU, S, t, logit.PI = TRUE){
  K <- length(PI)
  N <- dim(x)[1]
  p <- dim(x)[2]

  mp <- p * (p + 1) / 2
  M <- K - 1 + K * p + K * p * (p + 1) / 2
  G <- Gmat(p)
  Id <- diag(p)

  invS <- apply(S, 3, solve)
  dim(invS) <- c(p, p, K)

  SS <- lapply(1:N, function(i){
    Sbuf.1 <- NULL
    if(K != 1){
      if(logit.PI){    # partial wrt the multivariate logit parameter
        Sbuf.1 <- t[i, -K] - PI[-K]
      } else{
        tmp <- t[i, ] / PI
        Sbuf.1 <- tmp[1:(K-1)] - tmp[K]
      }
    }

    Sbuf.2 <- list()
    Sbuf.3 <- list()
    for(j in 1:K){
      tmp <- x[i,] - MU[j,]
      dim(tmp) <- c(p, 1)
      tmp.1 <- t[i,j] * invS[,,j]
      Sbuf.2[[j]] <- tmp.1 %*% tmp
      Smat <-  tmp.1 %*% (tmp %*% t(tmp) %*% invS[,,j] - Id) / 2
      Sbuf.3[[j]] <- as.vector(Smat) %*% G
    }

    c(Sbuf.1, do.call("c", Sbuf.2), do.call("c", Sbuf.3))
  })

  do.call("cbind", SS)
} # End of partial.q().


### Back compartiable with Volodymyr's code.
Iy <- function(x, PI, MU, S, t){
  SS <- partial.q(x, PI, MU, S, t, logit.PI = FALSE)
  SS %*% t(SS)
} # End of Iy().

Iy2 <- function(x, PI.0, MU.0, S.0, t.0, PI.a, MU.a, S.a, t.a){
  SS <- rbind(partial.q(x, PI.0, MU.0, S.0, t.0, logit.PI = FALSE),
              partial.q(x, PI.a, MU.a, S.a, t.a, logit.PI = FALSE))
  SS %*% t(SS)
} # End of Iy2().


### Compute posterior.
postPI <- function(x, emobj){
  e.step(x, emobj)$Gamma
} # End of postPI().


### Generate dataset.
GenDataSet <- function(N, PI, MU, S){
  K <- dim(MU)[1]
  p <- dim(MU)[2]
  Nk <- drop(rmultinom(1, N, PI))

  id <- rep(1:K, Nk)
  Sample <- lapply(1:K, function(i){
    ### To avoid degeneration.
    if(Nk[i] > 0){
      ret <- mvrnorm(n = Nk[i], mu = MU[i,], Sigma = S[,,i])
    } else{
      ret <- NULL
    }
    ret
  })
  x <- do.call("rbind", Sample)

  list(x = x, id = id)
} # End of GenDataSet().

GenMixDataSet <- function(N, PI.0, MU.0, S.0, PI.a, MU.a, S.a, tau = 0.5){
  N.0 <- rbinom(1, N, c(tau, 1-tau))
  N.a <- N - N.0

  ret.0 <- GenDataSet(N.0, PI.0, MU.0, S.0)
  ret.a <- GenDataSet(N.a, PI.a, MU.a, S.a)

  list(x = rbind(ret.0$x, ret.a$x), id = c(ret.0$id, ret.a$id),
       hid = c(rep(0, N.0), rep(1, N.a)))
} # End of GenMixDataSet().


### Observed Information for DataSet.
w.2 <- function(x, emobj.0, emobj.a, tau = 0.5){
  f.0 <- sum(log(dmixmvn(x, emobj.0))) + log(tau)
  f.a <- sum(log(dmixmvn(x, emobj.a))) + log(1 - tau)
  ### Numerical unstable.
  # g <- exp(f.0) + exp(f.a)
  # c(exp(f.0) / g, exp(f.a) / g)
  pi.0 <- 1 / (exp(f.a - f.0) + 1)
  c(pi.0, 1 - pi.0)
} # End of w.2().


### Obtain parameters.
get.E.chi2 <- function(x, emobj.0, emobj.a, given = c("0", "a"), tau = 0.5,
    n.mc = 1000, verbose = TRUE){
  N <- nrow(x)
  p <- ncol(x)

  K.0 <- emobj.0$nclass
  PI.0 <- emobj.0$pi
  MU.0 <- emobj.0$Mu
  S.0 <- LTSigma2variance(emobj.0$LTSigma) 

  K.a <- emobj.a$nclass
  PI.a <- emobj.a$pi
  MU.a <- emobj.a$Mu
  S.a <- LTSigma2variance(emobj.a$LTSigma) 

  if(given[1] == "0"){
    PI <- PI.0
    MU <- MU.0
    S <- S.0
  } else if(given[1] == "a"){
    PI <- PI.a
    MU <- MU.a
    S <- S.a
  } else{
    stop("given should be '0' or 'a'.")
  }

  ### Obtain nabla logL via Monte Carlo.
  x.new <- GenDataSet(n.mc, PI, MU, S)$x
  t.0 <- postPI(x.new, emobj.0)
  t.a <- postPI(x.new, emobj.a)
  pl.0 <- partial.q(x.new, PI.0, MU.0, S.0, t.0, logit.PI = FALSE)
  pl.a <- partial.q(x.new, PI.a, MU.a, S.a, t.a, logit.PI = FALSE)

  ### Expected degrees of freedom.
  nl <- rbind(pl.0, pl.a)
  mu <- rowMeans(nl)
  nl <- nl - mu
  J <- nl %*% t(nl) / n.mc
  par.df <- eigen(J, TRUE, only.values = TRUE)$values
  par.df <- par.df[par.df > 0]
  par.df <- sum((par.df / par.df[1] > 1e-10) |
                (cumsum(par.df) / sum(par.df) < 0.90))

  ### Expected noncenteriality
  M.0 <- nrow(pl.0)
  M.a <- nrow(pl.a)
  if(given == "0"){
    id <- M.0 + (1:M.a)
  } else{
    id <- 1:M.0
  }
  J.nc <- matrix(J[id, id], nrow = length(id))
  nu <- matrix(mu[id], nrow = length(id))
  ### Numerical unstable.
  # par.nc <- t(nu) %*% solve(J.nc) %*% nu
  tmp <- eigen(J.nc, TRUE)
  tmp.nc <- tmp$values[tmp$values > 0]
  tmp.nc <- sum((tmp.nc / tmp.nc[1] > 1e-8) |
                (cumsum(tmp.nc) / sum(tmp.nc) < 0.90))
  nu <- t(nu) %*% tmp$vectors[, 1:tmp.nc]
  par.nc <- nu %*% diag(1 / tmp$values[1:tmp.nc], tmp.nc, tmp.nc) %*% t(nu)

  ### For returns.
  ret <- c(par.df, par.nc)
  if(verbose){
    cat("K.0=", K.0, ", M.0=", M.0, " v.s. K.a=", K.a, ", M.a=", M.a,
        " | given=", given, " : df=", ret[1], ", nc=", ret[2],
        ".\n", sep = "")
  }
  ret
} # End of get.E.chi2().

get.E.delta <- function(x, emobj.0, emobj.a, tau = 0.5, n.mc = 1000){
  N <- nrow(x)

  PI.0 <- emobj.0$pi
  MU.0 <- emobj.0$Mu
  S.0 <- LTSigma2variance(emobj.0$LTSigma) 

  PI.a <- emobj.a$pi
  MU.a <- emobj.a$Mu
  S.a <- LTSigma2variance(emobj.a$LTSigma) 

  ### This is inaccurate and incorrect!!!
  # E.delta <- lapply(1:n.mc, function(i){
  #   x.new <- GenMixDataSet(N, PI.0, MU.0, S.0, PI.a, MU.a, S.a, tau = tau)$x
  #   logL(x.new, emobj.a) - logL(x.new, emobj.0)
  # })

  n.mc.0 <- rbinom(1, n.mc, c(tau, 1-tau))
  n.mc.a <- n.mc - n.mc.0

  E.delta <- NULL
  if(n.mc.0 > 0){
    tmp <- lapply(1:n.mc.0, function(i){
      x.new <- GenDataSet(N, PI.0, MU.0, S.0)$x
      logL(x.new, emobj.a) - logL(x.new, emobj.0)
    })
    E.delta <- c(E.delta, tmp)
  }
  if(n.mc.a > 0){
    tmp <- lapply(1:n.mc.a, function(i){
      x.new <- GenDataSet(N, PI.a, MU.a, S.a)$x
      logL(x.new, emobj.a) - logL(x.new, emobj.0)
    })
    E.delta <- c(E.delta, tmp)
  }

  do.call("sum", E.delta) / n.mc
} # End of get.E.delta().


pchisq.my <- function(q, df, ncp = 0, lower.tail = TRUE){
  if(ncp == Inf){
    if(lower.tail){
      ret <- 0
    } else{
      ret <- 1
    }
  } else{
    ret <- pchisq(q, df, ncp, lower.tail)
  }
  ret
} # End of pchisq.my().
