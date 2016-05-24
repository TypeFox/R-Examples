rq.fit.panel <- function(X,y,s,w=c(.25,.5,.25),taus=(1:3)/4,lambda = 1){
# prototype function for panel data fitting of QR models
# the matrix X is assumed to contain an intercept
# the vector s is a strata indicator assumed (so far) to be a one-way layout
# NB:  
# 0.  This is an altered version from that originally posted -- the definition
#     of the rhs vector now incorporates a factor of 1/2 for the penalty.
# 1.  The value of the shrinkage parameter lambda is an open research problem in
#   the simplest homogneous settings it should be the ratio of the scale parameters
#   of the fixed effects and the idiocyncratic errors
# 2.  On return the coefficient vector has m*p + n elements where m is the number
# quantiles being estimated, p is the number of colums of X, and n is the
# number of distinct values of s.  The first m*p coefficients are the 
# slope estimates, and the last n are the "fixed effects"
# 3.  Like all shrinkage (regularization) estimators, asymptotic inference is somewhat
# problematic... so the bootstrap is the natural first resort.
  require(SparseM)
  require(quantreg)
  K <- length(w)
  if(K != length(taus))
    stop("length of w and taus must match")
  X <- as.matrix(X)
  p <- ncol(X)
  n <- length(levels(as.factor(s)))
  N <- length(y)
  if(N != length(s) || N != nrow(X))
    stop("dimensions of y,X,s must match")
  Z <- as.matrix.csr(model.matrix(~as.factor(s)-1))
  Fidelity <- cbind(as(w,"matrix.diag.csr") %x% X,w %x% Z)
  Penalty <- cbind(as.matrix.csr(0,n,K*p),lambda*as(n,"matrix.diag.csr"))
  D <- rbind(Fidelity,Penalty)
  y <- c(w %x% y,rep(0,n))
  a <- c((w*(1-taus)) %x% (t(X)%*%rep(1,N)),
    sum(w*(1-taus)) * (t(Z) %*% rep(1,N)) + lambda * rep(1/2,n))
  rq.fit.sfn(D,y,rhs=a)
}
