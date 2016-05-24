# This function performs cross-validation for Grace
# Author: Sen Zhao
# Email: sendavid7@gmail.com

cvGrace <- function(X, Y, L, lambda.L, lambda.1, lambda.2, K = 10){
  lambda.1 <- sort(lambda.1, decreasing = TRUE)
  lambda.L <- sort(lambda.L, decreasing = TRUE)
  lambda.2 <- sort(lambda.2, decreasing = TRUE)
  if(length(lambda.1) == 1){
    lambda.1 <- rep(lambda.1, 2)
  }
  p <- ncol(X)
  n <- nrow(X)
  lam1 <- matrix(nrow = length(lambda.L), ncol = length(lambda.2))
  ERR <- matrix(0, nrow = length(lambda.L), ncol = length(lambda.2))
  for(iL in 1:length(lambda.L)){
    lL <- lambda.L[iL]
    for(i2 in 1:length(lambda.2)){
      l2 <- lambda.2[i2]
      Lnew <- lL * L + l2 * diag(p)
      S <- eigen(Lnew)$vectors %*% sqrt(diag(eigen(Lnew)$values))
      l2star <- 1
      l1star <- lambda.1
      Xstar <- rbind(X, sqrt(l2star) * t(S)) / sqrt(1 + l2star)
      Ystar <- c(Y, rep(0, p))
      gammastar <- l1star / sqrt(1 + l2star) / 2 / (n + p)
      cvres <- cv.glmnet(X, Y, lambda = gammastar, intercept = FALSE, standardize = FALSE, nfolds = K)
      lam1[iL, i2] <- which.min(abs(gammastar - cvres$lambda.min))
      ERR[iL, i2] <- cvres$cvm[lam1[iL, i2]]
    }
  }
  idx <- which(ERR == min(ERR), arr.ind = TRUE)
  resl1 <- lambda.1[lam1[idx]][1]
  reslL <- lambda.L[idx[1, 1]]
  resl2 <- lambda.2[idx[1, 2]]
  return(c(reslL, resl1, resl2))
}