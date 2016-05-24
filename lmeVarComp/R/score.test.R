score.test <-
function(Y, X, Z, Sigma, m0)
{
  require(varComp)

  m1 <- length(Z)
  K <- lapply(seq_len(m1), function(j) 
    tcrossprod(Z[[j]] %*% Sigma[[j]], Z[[j]]))
  fit <- varComp::varComp.fit(Y, X, K[seq_len(m0)])
  test <- varComp::varComp.test(fit, additional.varcov = K[(m0 + 1L) : m1])
  result <- test[[1L]][[1L]][[1L]]
  c(stat.obs = result$statistic, p.value = result$p.value)
}
