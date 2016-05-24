bayes.regress <- function(formula, data, sample = 1000, weights =
                          NULL, truncate=FALSE){
  fml <- as.formula(formula)
  D <- model.frame(fml, data = data)
  Y <- as.matrix(model.response(D))
  X <- model.matrix(formula, data = D)
  if (!is.null(weights)) {
    weights <- weights/sum(weights)
    w <- weights^(-1/2)
    X <- w * X
    Y <- w * Y
  }
  n <- nrow(X)
  nu <- n - ncol(X)
  s2 <- as.numeric((t(Y) %*% Y - t(Y) %*% X %*% solve(t(X) %*% X) %*% t(X) %*% Y)) / nu
  phi <- rchisq(sample, df = nu)
  sigma2 <- (nu * s2) / phi
  beta.hat <- solve(t(X) %*% X) %*% t(X) %*% Y
  tXX <- function(sigma2, X) sigma2 * solve(t(X) %*% X)
  Sigma <- lapply(sigma2, tXX, X)

  truncMVRnorm <- function(sig, mean, n){
    repeat{
      draw <- mvrnorm(n=n, mu=mean, Sigma=sig)
      if(min(draw)>=0 & max(draw)<=1){
       
        return(draw)
      }
    }
  }
  
  if(truncate==TRUE){
    beta <- t(sapply(Sigma, truncMVRnorm, n = 1, mean = beta.hat))
  }else{
    beta <- t(sapply(Sigma, mvrnorm, n = 1, mu = beta.hat))
  }
  colnames(beta) <- colnames(X)
  beta
}
