logLaplace <-
function (formula, data, tools) {

  formula <- as.formula(formula)
  X <- model.matrix(formula, data)
  y <- data$freq
  fit <- glm (formula, data, family = poisson())
  W <- diag (fit$fitted.values)
  maxloglik <- sum(y  * log(fit$fitted.values))
  X <- X[,-1]
  fisherMatrix <- t(X) %*% W %*% X
  J <- dim(X)[2]
  d <- det(fisherMatrix)
  
  if (d == 0) {
    model <- findModel(formula, tools)
    generators <- findGenerators (model, tools)
    stop(paste ("Hessian is numerically singular for model: \n", prettyHierFormula(generators, tools), ".\n Hierarchical log-linear model search can not be completed.\n", "Try increasing alpha.", sep = ""))
  }

  value <- maxloglik + 0.5 * J * log (2*pi) - 0.5 * log(d) 
  
  return(value)   

}
