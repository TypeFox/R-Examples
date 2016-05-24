pseudo.rlr.test <-
function(Y, X, Z, Sigma, m0, nsim = 5000L, seed = 130623L)
{
  require(varComp)
  require(RLRsim)

  j0 <- seq_len(m0)
  j1 <- (m0 + 1L) : length(Z)
  m1 <- m0 + 1L
  
  Z[[m1]] <- do.call(cbind, Z[j1])
  Sigma[[m1]] <- do.call(block.diag, Sigma[j1])
  
  K <- lapply(seq_len(m1), function(j) 
    tcrossprod(Z[[j]] %*% Sigma[[j]], Z[[j]]))
  fit1 <- varComp::varComp.fit(Y, X, K)
  fit0 <- varComp::varComp.fit(Y, X, K[j0])
  
  stat.obs <- 2 * (fit1$PREML - fit0$PREML)
  if (stat.obs < qchisq(0.2, 1)) {
    p.value <- 1
  } else if (stat.obs > qchisq(0.998, 1)) {
    p.value <- 0
  } else if (nsim >= 1L) {
    # avoids chol(Sigma[[m1]])
    decomp <- eigen(Sigma[[m1]], TRUE)
    v <- decomp$values
    v[v < 0] <- 0
    R <- decomp$vectors %*% (sqrt(v) * t(decomp$vectors))
    stat.sim <- RLRsim::RLRTSim(X, Z[[m1]], qr(X), R, nsim = nsim, seed = seed)
    p.value <- mean(stat.sim >= stat.obs)
  } else {
    p.value <- NA_real_
  }
  
  c(stat.obs = stat.obs, p.value = p.value)
}
