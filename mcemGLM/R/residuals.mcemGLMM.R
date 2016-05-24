residuals.mcemGLMM <- function(object, type = c("deviance", "pearson"), ...) {
  kP <- ncol(object$x)
  coef0 <- tail(object$mcemEST, 1)[1:kP]
  u0 <- colMeans(object$randeff)
  lin0 <- object$x %*% coef0 + object$z %*% u0
  if (object$call$family == "bernoulli") {
    prob0 <- exp(lin0) / (1 + exp(lin0))
    if (type[1] == "pearson") {
      res0 <- (object$y - prob0) / sqrt(prob0 * (1 - prob0))
    }
    if (type[1] == "deviance") {
      res0 <- ifelse(object$y, 1, -1) * sqrt(-2 * (object$y * log(prob0) + (1 - object$y) * log(1 - prob0)))
    }
  }
  
  if (object$call$family == "poisson") {
    mu0 <- exp(lin0)
    if (type[1] == "pearson") {
      res0 <- (object$y - mu0) / sqrt(mu0)
    }
    if (type[1] == "deviance") {
      # print(mu0)
      res0 <- sign(object$y - mu0) * sqrt(2 * ifelse(object$y > 0, object$y * log(object$y / mu0), 0) - 2 * (object$y - mu0))
    }
  }
  
  if (object$call$family == "negbinom") {
    mu0 <- exp(lin0)
    a0 <- tail(object$mcemEST, 1)[kP + 1]
    if (type[1] == "pearson") {
      res0 <- (object$y - mu0) / sqrt(mu0 * (1 + 1/a0))
    }
    if (type[1] == "deviance") {
      res0 <- sign(object$y - mu0) * sqrt(2 * (ifelse(object$y > 0, object$y * log(object$y / mu0), 0)) - 2 * (object$y + a0) * log((object$y + a0)/(mu0 + a0)))
    }
  }
  
  if (object$call$family == "gamma") {
    mu0 <- exp(lin0)
    a0 <- tail(object$mcemEST, 1)[kP + 1]
    if (type[1] == "pearson") {
      res0 <- (object$y - mu0) / (mu0 / sqrt(a0))
    }
    if (type[1] == "deviance") {
      res0 <- sign(object$y - mu0) * sqrt(-2 * a0 * (log(object$y/mu0) - (object$y - mu0)/mu0))
    }
  }
  return(as.vector(res0))
}