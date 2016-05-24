fit.enet.fixed.lambda <-
function (y, X, lambda, alpha, mu){
  
  ### first sort out the standardisations of y and X, based on the specification of mu
  if (is.null(mu)){ # mean center y and columns of X
    y.tilde = y - mean(y)
    X.tilde = apply (X, 2, function(col){col-mean(col)})
  }else{ # just subtract mu from y
    y.tilde = y - mu
    X.tilde = X
  }
  
  ## first get the sequence at which glmnet would fit the model
  max.lambda = max (abs(t(X.tilde)%*%y.tilde))
  #print (max.lambda)
  lambda.seq = sort(c(lambda, exp(seq (from=log(max.lambda), to=log(max.lambda) - 4*log(10), length=100))), decreasing=TRUE)
  which.ind = which (lambda.seq==lambda)
  lambda.seq = lambda.seq/length(y) # scale befre call to glmnet
  
  glmnet.obj = glmnet (y=y.tilde, x=X.tilde, standardize=FALSE, intercept=FALSE, alpha=alpha, lambda=lambda.seq, thresh=10^-15, maxit=10^7)
  beta = glmnet.obj$beta[,which.ind]
  
  list(beta=beta, lambda=lambda, alpha=alpha, X=X.tilde, y=y.tilde, mu=mu)
}
