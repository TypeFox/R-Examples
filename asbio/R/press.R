press <- function(lm){
X <- model.matrix(lm)
H <- X%*%(solve((t(X)) %*% (X))) %*% t(X)
e.i <- resid(lm)
sum((e.i/(1-diag(H)))^2)
}

