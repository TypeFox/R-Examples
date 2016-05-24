logit.weight.test <-
function(y,X,G,Z,maf,weight.beta=c(1,25),method="liu")
{
  library(CompQuadForm)

  y <- as.vector(y)
  G <- as.matrix(G)
  X <- as.matrix(X)
  Z <- as.matrix(Z)

  GZ <- G %*% Z
  M <- cbind( X, GZ )

  ## fit null model

  ## M_0: logit P(Y=1)  = X\alpha 
  fit.0 <- glm( formula = y ~ X-1, family = binomial( link = logit ) )
  mu.0 <- fit.0$fitted.value
  d.0 <- mu.0*(1-mu.0)
  res.0 <- y - mu.0

  ## M_{0a}: logit P(Y=1) = X\alpha + GZ\pi 
  fit.0a <- glm( formula = y ~ -1 + X + GZ, family = binomial( link = logit ) )
  mu.0a <- fit.0a$fitted.value
  d.0a <- mu.0a*(1-mu.0a)
  res.0a <- y - mu.0a

  ## define useful matrices
  
  n <- dim(X)[1]
  I <- diag(1,n) 

  D.0 <- diag(d.0)
  D.0a <- diag(d.0a)

  ## ( X^T %*% D.0 %*% X )^{-1}
  tXD0X <- t(X) %*% D.0 %*% X
  inv.tXD0X <- solve( tXD0X )

  ## ( M^T %*% D.0a %*% M )^{-1}
  tMD0aM <- t(M) %*% D.0a %*% M
  inv.tMD0aM <- solve( tMD0aM )

  #P1 <- I - D.0 %*% X %*% inv.tXD0X %*% t(X) 
  #P2 <- I - D.0a %*% M %*% inv.tMD0aM %*% t(M)

  P01 = D.0 - (d.0 %o% d.0) * ( X %*% (inv.tXD0X) %*% t(X) )
  P02 = D.0a - (d.0a %o% d.0a) * ( M %*% (inv.tMD0aM) %*% t(M) ) 

  ## find square root of matrices P01 and P02
  #P01.half <- root.mat( P01 )
  #P02.half <- root.mat( P02 )


  ## calculate test statistics

  W <- diag( dbeta(maf,weight.beta[1],weight.beta[2])^2 )
  
  S.tau <-  0.5 * t(res.0a) %*% G %*% W %*% t(G) %*% res.0a

  inv.I.pi <- solve( t(GZ) %*% P01 %*% GZ )

  S.pi <- t(res.0) %*% GZ %*% inv.I.pi %*% t(GZ) %*% res.0

  ## null distributions and p-values

  ## for S.pi
  p.value.S.pi <- 1-pchisq(S.pi,df=dim(Z)[2])

  ## for S.tau
  #P2.half.G <- P02.half %*% G %*% sqrt(W)
  #Mat <- 0.5 * t( P2.half.G ) %*% P2.half.G
  Mat <- 0.5* sqrt(W) %*% t(G) %*% P02 %*% G %*% sqrt(W)
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

