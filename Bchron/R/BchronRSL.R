BchronRSL <-
function(BchronologyRun,RSLmean,RSLsd,degree=1,iterations=10000,burn=2000,thin=8) {
  if(degree>5) stop('Degree not supported')
  remaining=(iterations-burn)/thin
  betaStore = matrix(NA,ncol=degree+1,nrow=remaining)
  y = matrix(RSLmean,ncol=1)
  Q = diag(1/((RSLsd)^2))
  N = nrow(y)
  chrons = BchronologyRun$thetaPredict/1000
  whichrows = sample(1:nrow(chrons),iterations,replace=TRUE)
  degmat = matrix(rep(0:(degree),ncol(chrons)*N),nrow=N,ncol=degree+1,byrow=TRUE)
  const = mean(as.matrix(chrons))
  pb = utils::txtProgressBar(min = 1, max = iterations, style = 3,width=60,title='Running BchronRSL')
  for(i in 1:iterations) {
    utils::setTxtProgressBar(pb, i)
    if(i%%20==0 | i==1) {
      currchron = chrons[whichrows[i],] - const
      X = matrix(rep(as.numeric(currchron),degree+1),ncol=degree+1)
      X = X^degmat
    }
    # Sample a beta
    if(i%%thin==0 & i>burn) {
      betaStore[(i-burn)/thin,] = matrix(MASS::mvrnorm(1,solve(t(X)%*%Q%*%X,t(X)%*%Q%*%y),solve(t(X)%*%Q%*%X)),ncol=1)
    }
  }
    
  out = list(BchronologyRun=BchronologyRun,samples=betaStore,degree=degree,RSLmean=RSLmean,RSLsd=RSLsd,const=const)
  class(out) = 'BchronRSLRun'
  return(out)
  
}
