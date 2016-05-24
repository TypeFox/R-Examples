Kappa <-
function(beta, X, Y, lambda_kappa){
  n <- nrow(Y)  
  SigmaR <- 1/n * t(Y - X %*% beta) %*% (Y - X %*% beta)
  res <- glasso(SigmaR, lambda_kappa, penalize.diagonal = FALSE)
  return(as.matrix(forceSymmetric(res$wi)))
}
