error.nsc <-function(object) {
###Computes the roc curve for a nsc model
  yhat <- object$yhat
  y <- object$y
  ny <- table(y)
  errors <- matrix(0, length(object$threshold), length(ny))
  Y <- data.matrix(yhat) != unclass(y)
  yind <- model.matrix( ~ factor(y) - 1, data = list(y = y))
  errors <- t(t(yind) %*% Y)
  apply(errors, 2, mean)
}
