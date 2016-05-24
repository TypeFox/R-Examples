`Stem.Simulation` <-
function (StemModel){

  require(MASS)

  n=StemModel$data$n
  d=StemModel$data$d
  r=StemModel$data$r
  p=StemModel$skeleton$p
  
  z = matrix(NA,nrow = d,ncol = n)
  y = matrix(NA,nrow = p,ncol = n)

  phi 			= StemModel$skeleton$phi
  phi$beta		= matrix(phi$beta,r,1)
  phi$logb 		= log(phi$sigma2eps/phi$sigma2omega)
  phi$logtheta	= log(phi$theta)

  dist     = as.matrix(dist(StemModel$data$coordinates,diag=TRUE)) #distance matrix
  covariates 	= StemModel$data$covariates
  covariates  = changedimension_covariates(covariates,d=d,r=r,n=n)

  ####################
  ###Matrix definitions
  ####################
  Fmat = StemModel$skeleton$K
  Gmat = phi$G                     

  #cov.spaz = phi_real$sigma2omega * exp(-phi_real$theta * dist)
  #Vmat     = diag(phi_real$sigma2eps,d) + cov.spaz
  
  Vmat=phi$sigma2omega * Sigmastar.exp(d=d,logb=phi$logb,logtheta=phi$logtheta,dist=dist)
  
  Wmat = phi$Sigmaeta                           

  m0   = phi$m0                                 
  C0   = phi$C0		                     

  #####################
  y0 = matrix(mvrnorm(n=1, mu=m0, Sigma=C0),nrow=p, ncol=1)

  #t=1
  y[,1] = Gmat %*% y0 + mvrnorm(n=1, mu=matrix(0,nrow=p,ncol=1), Sigma=Wmat)
  z[,1] = covariates[,,1] %*% phi$beta + Fmat %*% y[,1] + mvrnorm(n=1, mu=matrix(0,nrow=d,ncol=1),Sigma = Vmat)

  for (t in 2:n) {
      y[,t] = Gmat %*% y[,t-1] + mvrnorm(n=1, mu=matrix(0,nrow=p,ncol=1), Sigma=Wmat)
      z[,t] = covariates[,,t] %*% phi$beta + Fmat %*% y[,t]   + mvrnorm(n=1, mu=matrix(0,nrow=d,ncol=1),Sigma=Vmat)
  }
  

  rownames(z) = rownames(StemModel$data$coordinates)
  return(z=t(z))
}

