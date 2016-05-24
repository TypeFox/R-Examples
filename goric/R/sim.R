sim <- function(object, n.sims){
  UseMethod("sim")
}

sim.orlm <- function(object, n.sims){
  object <- object
  unr <- lm(object$y ~ object$X-1)
  ### code taken from package arm
  summ <- summary(unr)
  coef <- summ$coef[,1:2,drop=FALSE]
  dimnames(coef)[[2]] <- c("coef.est","coef.sd")
  sigma.hat <- summ$sigma
  beta.hat <- coef[,1,drop = FALSE]
  V.beta <- summ$cov.unscaled
  n <- summ$df[1] + summ$df[2]
  k <- summ$df[1]
  sigma <- sigma.hat*sqrt((n-k)/rchisq(n.sims,n-k))
  # new simulated dataset
  smat <- matrix(rep(object$X %*% coefficients(unr), n.sims), ncol=n.sims) + sapply(sigma, function(sigma) rnorm(nrow(object$X), 0, sigma))
  # new set of constraint parameters
  beta <- rbind(apply(smat, 2, function(sy) coefficients(orlm(sy ~ object$X-1, constr=object$constr, rhs=object$rhs, nec=object$nec))))
  rownames(beta) <- colnames(object$X)
  return(list(coef=beta, sigma=sigma))
}

sim.orgls <- function(object, n.sims){
  object <- object
  unr <- lm(object$y ~ object$X-1)
  ### code taken from package arm
  summ <- summary(unr)
  coef <- summ$coef[,1:2,drop=FALSE]
  dimnames(coef)[[2]] <- c("coef.est","coef.sd")
  sigma.hat <- summ$sigma
  beta.hat <- coef[,1,drop = FALSE]
  V.beta <- summ$cov.unscaled
  n <- summ$df[1] + summ$df[2]
  k <- summ$df[1]
  sigma <- sigma.hat*sqrt((n-k)/rchisq(n.sims,n-k))
  # new simulated dataset
  smat <- matrix(rep(object$X %*% coefficients(unr), n.sims), ncol=n.sims) + sapply(sigma, function(sigma) rnorm(nrow(object$X), 0, sigma))
  # new set of constraint parameters
  beta <- rbind(apply(smat, 2, function(sy) coefficients(orlm(sy ~ object$X-1, constr=object$constr, rhs=object$rhs, nec=object$nec))))
  rownames(beta) <- colnames(object$X)
  return(list(coef=beta, sigma=sigma))
}
