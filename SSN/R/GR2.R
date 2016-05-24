GR2 <-
function(object)
{
  if(class(object) != "glmssn") return("Not a glmssn object")
  W <- object$sampinfo$X
  Vi <- object$estimates$Vi
  z <- object$sampinfo$z
  betahat <- object$estimates$betahat
  muhat <- sum(Vi %*% z)/sum(Vi)
  ones <- matrix(1, ncol = 1, nrow = length(z))
  1 - t(z - W %*% betahat) %*% Vi %*% (z - W %*% betahat)/
    t(z - ones %*% muhat) %*% Vi %*% (z - ones %*% muhat)
}

