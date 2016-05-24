# This function calculates Grace coefficient estimates
# Author: Sen Zhao
# Email: sendavid7@gmail.com

grace <- function(Y, X, L, lambda.L, lambda.1 = 0, lambda.2 = 0, normalize.L = FALSE, K = 10){
  ori.Y <- Y
  ori.X <- X
  if(length(Y) != nrow(X)){
    stop("Error: Dimensions of X and Y must match.")
  }
  if(!isSymmetric(L)){
    stop("Error: L must be a symmetric matrix.")
  }
  if(ncol(X) != ncol(L)){
    stop("Error: Dimensions of X and L must match.")
  }
  if(min(lambda.L) < 0 | min(lambda.2) < 0){
    stop("Error: Grace tuning parameters must be non-negative.")
  }
  if(min(lambda.L) <= 0 & min(lambda.2) <= 0){
    stop("Error: At least one of the grace tuning parameters must be positive.")
  }
  
  Y <- Y - mean(Y)  # Center Y
  n <- nrow(X)
  p <- ncol(X)
  scale.fac <- attr(scale(X), "scaled:scale")
  X <- scale(X)     # Standardize X
  
  if(normalize.L){
    diag(L)[diag(L) == 0] <- 1
    L <- diag(1 / sqrt(diag(L))) %*% L %*% diag(1 / sqrt(diag(L)))  # Normalize L
  }
  
  emin <- min(eigen(lambda.L * L + lambda.2 * diag(p))$values)  # Minimum eigenvalue of the penalty weight matrix
  if(emin < 1e-5){
    lambda.2 <- 1e-5 - emin
    warning(paste("Warning: The penalty matrix (lambda.L * L + lambda.2 * I) is not positive definite. lambda.2 is adjusted to ", round(lambda.2, 5), " to make it positive definite.", sep = ""))
  }
  

  # If the more than one tuning parameter is provided, perform K-fold cross-validation  
  if((length(lambda.L) > 1) | (length(lambda.1) > 1) | (length(lambda.2) > 1)){
    tun <- cvGrace(X, Y, L, lambda.L, lambda.1, lambda.2, K)
    lambda.L <- tun[1]
    lambda.1 <- tun[2]
    lambda.2 <- tun[3]  
  }
  
  Lnew <- lambda.L * L + lambda.2 * diag(p)
  S <- eigen(Lnew)$vectors %*% sqrt(diag(eigen(Lnew)$values))
  l2star <- 1
  l1star <- lambda.1
  Xstar <- rbind(X, sqrt(l2star) * t(S)) / sqrt(1 + l2star)
  Ystar <- c(Y, rep(0, p))
  gammastar <- l1star / sqrt(1 + l2star) / 2 / (n + p)
  betahatstar <- glmnet(Xstar, Ystar, lambda = gammastar, intercept = FALSE, standardize = FALSE, thresh = 1e-11)$beta[, 1]
  betahat <- betahatstar / sqrt(1 + l2star)
  
  truebetahat <- betahat / scale.fac  # Scale back coefficient estimate
  truealphahat <- mean(ori.Y - ori.X %*% truebetahat)
  return(list(intercept = truealphahat, beta = truebetahat))
}
