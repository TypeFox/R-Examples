minque0 <-
function(Y, D, U, m0, m1)
{
  Ytilde2 <- as.double(crossprod(U, Y)) ^ 2
  
  # computes scaling factors for columns of D
  Sa <- crossprod(D)
  Sb <- colSums(D)
  S <- rbind(cbind(Sa, Sb), c(Sb, nrow(D)))
  Q <- c(crossprod(D, Ytilde2), sum(Ytilde2))
  
  k <- m1 + 1L
  tau <- mnls(S, Q)
  scaling <- pmax(tau[-k] / tau[k], 
    apply(D, 2L, function(v) 1 / quantile(v[v > 0], 0.99)))
  for (k in seq_len(m1)) {
    D[, k] <- D[, k] * scaling[k]
  }
  
  # computes initial estimates of tau1 and tau0
  Sa <- crossprod(D)
  Sb <- colSums(D)
  S <- rbind(cbind(Sa, Sb), c(Sb, nrow(D)))
  Q <- c(crossprod(D, Ytilde2), sum(Ytilde2))
  
  k <- m1 + 1L
  S1 <- mnls(S, diag(k))
  S1 <- (S1 + t(S1)) / 2
  tau1 <- c(S1 %*% Q)
  tau1 <- tau1[-k] / tau1[k]
  tau1[tau1 < 0] <- 0
  
  k <- m0 + 1L
  i0 <- c(seq_len(m0), m1 + 1L)
  S0 <- mnls(S[i0, i0, drop = FALSE], diag(k))
  S0 <- (S0 + t(S0)) / 2
  tau0 <- c(S0 %*% Q[i0])
  tau0 <- tau0[-k] / tau0[k]
  tau0[tau0 < 0] <- 0
  
  list(Ytilde2 = Ytilde2, D = D, scaling = scaling,
    S1 = S1, S0 = S0, tau1 = tau1, tau0 = tau0)
}
