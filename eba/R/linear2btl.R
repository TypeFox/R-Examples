linear2btl <- function(object, order=FALSE){
  # Takes a lm or glm object and computes BTL parameters
  #   and the covariance matrix using the delta method.
  # Assumes a design matrix as returned by pcX() with the reference
  #   category being the first level.
  # Discards the order effect when computing the covariance matrix.
  # 
  # Author: Florian Wickelmaier <wickelmaier@web.de>
  # Last mod: 2007/Oct/04, FW

  if(order){
    beta <- object$coef[-length(object$coef)]
    cov.beta <- summary(object)$cov.un[1:length(beta), 1:length(beta)]
  }else{
    beta <- object$coef
    cov.beta <- summary(object)$cov.un
  }

  estimate <- exp(c(0, beta)) / sum(exp(c(0, beta)))  # BTL parameters
  h        <- estimate[1] * estimate[-1]
  H1       <- estimate[-1] %*% t(estimate[-1]) - diag(estimate[-1])
  H        <- rbind(h, H1)
  cov.btl  <- H %*% cov.beta %*% t(H)
  out <- cbind(estimate, se=sqrt(diag(cov.btl)))
  ans <- list(btl.parameters=out, cova=cov.btl, linear.coefs=beta)
  ans
}
