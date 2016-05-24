linear.test <-
function(y,X,G,Z,method="liu")
## y is a vector of phenotype; X is the design matrix for 1st level covariates, including intercept;
## G is genotype matrix; Z is the design matrix for the 2nd level covariates, including intercept;
{
  library(CompQuadForm)

  y <- as.vector(y)
  G <- as.matrix(G)
  X <- as.matrix(X)
  Z <- as.matrix(Z)

  GZ <- G %*% Z
  M <- cbind( X, GZ )

  ## ( X^T %*% X )^{-1}
  tXX <- t(X) %*% X
  inv.tXX <- solve(tXX)

  ## ( M^T %*% M )^{-1}
  tMM <- t(M) %*% M
  inv.tMM <- solve(tMM)

  ## fit null model

  ## M_0: Y  = X\alpha + \epsilon
  fit.0 = lm( y ~ X -1 ) 
  tilde.sigma2 = summary(fit.0)$sigma^2
  res.0=fit.0$resid

  ## M_{0a}: Y = X\alpha + GZ\pi + \epsilon
  fit.0a = lm( y ~ M -1 )
  hat.sigma2 = summary(fit.0a)$sigma^2
  res.0a=fit.0a$resid 


  ## define useful matrices
  
  n <- dim(X)[1]
  I <- diag(1,n) 

  P2 <- I - X %*% inv.tXX %*% t(X) 
  P1 <- I - M %*% inv.tMM %*% t(M)

  ## calculate test statistics
  
  S.tau <-  t(res.0a) %*% G %*% t(G) %*% res.0a

  inv.I.pi <- solve( t(GZ) %*% P2 %*% GZ )

  S.pi <- t(res.0) %*% GZ %*% inv.I.pi %*% t(GZ) %*% res.0
  S.pi <- S.pi / tilde.sigma2

  ## null distributions and p-values

  ## for S.pi
  p.value.S.pi <- 1-pchisq(S.pi,df=dim(Z)[2])

  ## for U.tau.hatpi
  P1.G <- P1 %*% G
  Mat <- ( hat.sigma2 ) * t( P1.G ) %*% P1.G
  eigen.value <- eigen( Mat, symmetric=TRUE )$values
  lambda <- eigen.value

  if( method == "davies" )    
      {  p.value.S.tau <- davies( S.tau , lambda )$Qq }
  if( method == "liu" )    
      {  p.value.S.tau <- liu( S.tau , lambda ) }

  ## overall p-value: Fisher method
  q.fisher <- -2*( log( p.value.S.tau ) + log( p.value.S.pi ) )
  p.value.overall <- 1-pchisq( q.fisher, df=4 )

  ##### output #####
  out <-  list ( S.tau = S.tau, 
                 S.pi = S.pi,
                 p.value.S.pi = p.value.S.pi,
                 p.value.S.tau = p.value.S.tau,
                 p.value.overall = p.value.overall
               )
  return(out)
}

