SparseTSCGM_lambdas <-
function(X, Y, nlambda = 100){
  # From SparseTSCGM package:
  lambda.seq <- function(SS, SA,nlambda)
  {
    if (length(nlambda)==1) nlambda <- rep(nlambda,2)
    lambda.min.ratio=0.1
    d= dim(SS)[2]
    lambda.max1 = max(max(SS-diag(d)),-min(SS-diag(d)))
    lambda.min1 = lambda.min.ratio*lambda.max1
    lambda1 = exp(seq(log(lambda.max1), log(lambda.min1), length = nlambda[1]))
    lambda.min.ratio2=0.15
    lambda.max2 = max(max(SA),-min(SA))
    lambda.min2 = lambda.min.ratio2*lambda.max2
    lambda2 = exp(seq(log(lambda.max2), log(lambda.min2), length = nlambda[2]))
    return(list(lambda1=lambda1, lambda2=lambda2))
  }
  
  T <- dim(Y)[1]
  p <- dim(X)[2]
  n <- 1
  q <- dim(Y)[2]
  xtyi <- array(NA, c(p,q,n))
  xtxi <- array(NA, c(p,p,n))
  ytyi <- array(NA, c(q,q,n))
  
    XX <- X
    YY <- Y
    XX2 <- X^2
    YY2 <- Y^2
    xtyi <- crossprod(XX,YY)
    xtxi <- crossprod(XX)
    ytyi <- crossprod(YY)

  xty=apply(xtyi, c(1,2), sum)
  xtx=apply(xtxi, c(1,2), sum)
  yty=apply(ytyi, c(1,2), sum)
  xtxt=apply(xtxi, c(1,2), sum)/(n*T)
  xtx2=(n*T)*colMeans(apply(XX2, c(1,2), sum))
  yty2=(n*T)*colMeans(apply(YY2, c(1,2), sum))
  
  SX <- xtx/(n*T)
  mSX <- glasso(SX,0.05,penalize.diagonal=FALSE)
  
  SX <- xtx/(n*T)
  mSX <- glasso(SX,0.05,penalize.diagonal=FALSE)
  SXi <- mSX$wi
  SS =(yty)/(n*T)
  SS = cov2cor(SS)
  SAs = xty/(n*T)
  SA = SAs %*% SXi
  
  lambda <-  lambda.seq(SS=SS,SA=SA, nlambda=nlambda)
  lam1 <- round(lambda$lambda1,3) 
  lam2 <- round(lambda$lambda2,3)
  lam2 <- round(lam2/max(lam2),3)    
  return(list(lambda_kappa = lam1, lambda_beta = lam2))
}
