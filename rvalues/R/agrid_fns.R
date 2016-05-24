## Functions in this file: rvalue.agrid.nn, rvalue.agrid.bb, rvalue.agrid.pg

###############################################################################
#######    r-values with conjugate/parametric priors  #########################
###############################################################################

###  normal/normal

rvalue.agrid.nn <- function( MLE.local, SE.local,hypers,alpha.grid,smooth)
{
  ## computes `R' values using normal posteriors and the `alpha grid' algorithm
  ## compared to Aug 2013 version, this one takes input on mu and tau2
  ## and it does an alpha grid that's regular on the log scale in (eps, 1/2]
  ##
  ###############################################################
  ## Input
  ##
  ## MLE.local  vector of local MLEs
  ## SE.local	matched vector of local standard errors
  ## mu           center of MLE.local (to be removed)
  ## tau2         prior variance of theta's 
  
  ################################################################
  ## input checks: TBD
  
  ################################################################
  ##
  ## of the conjugate normal(mu,tau2) hyperparams by approximate MLE; or input
  
  if( length(hypers)==2  )
  { mu <- hypers[1]; tau2 <- hypers[2] }
  if( hypers[1]=="estimate" )
  {
    mu <- mean(MLE.local)
    gethyps <- HyperEstimate(MLE.local,SE.local,family="gaussian")
    tau2 <- gethyps[2]  ## MLE of prior variance, pre-scaling
  }
  
  ################################################################
  ## standardize 
  
  xx <- (MLE.local - mu)/sqrt(tau2)  
  sig2 <- (SE.local^2)/tau2
  
  ################################################################
  ## get posterior means and SD's on this scale
  
  thetaPM <- xx/(1+sig2)
  thetaSD <- sqrt( sig2/(1+sig2) )
  
  ################################################################
  ## alpha grid algorithm: normal posteriors ; empirical marginals
  ##
  nunits <- length(xx)
  ngrid <- length(alpha.grid)
  
  thetaAlpha <- qnorm( alpha.grid, lower.tail=FALSE )
  
  dev <- (-1)*( outer( thetaPM, thetaAlpha, FUN="-" ) )
  z <- dev/thetaSD   ## nunits x ngrid matrix
  
  V <- pnorm( z, lower.tail=FALSE )  ## V[i,a] = P(theta_i>=theta_a|unit i data)
  ### Is it better to work on the log scale?
  
  tmp <- rvalueGuts( dat=cbind(xx,sig2), alpha.grid=alpha.grid, V=V, vfun=vfun.nn, 
                     hypers=c(0,1),smooth=smooth)
  rvals <- tmp$rvals
  
  bar <- data.frame( RValue=rvals, RV.rank=rank(rvals),
                     MLE.rank=rank(-MLE.local),
                     PM.rank=rank(-thetaPM),
                     MLE=MLE.local,
                     SE=SE.local, 
                     PostMean=mu+sqrt(tau2)*thetaPM,
                     PostSD=sqrt(tau2)*thetaSD,
                     PVal = pnorm( MLE.local/SE.local, lower.tail=FALSE ),
                     PVal.rank = rank( -MLE.local/SE.local ) )  
  
  ord <- order( rvals, -thetaPM )
  res <- bar[ord,]
  
  theta.quantiles <- qnorm(alpha.grid,mean=mu,sd=sqrt(tau2),lower.tail=FALSE)
  aux <- list( V=V,  alpha.grid=alpha.grid, Vmarginals=tmp$lamfun,Vmarginals.smooth=tmp$smoothlamfun,
               unsorted = bar,hypers=c(mu, tau2),theta.quantiles=theta.quantiles,
               prior="conjugate",family="gaussian",smooth=smooth)
  return( list( main=res, aux=aux, rvalues=rvals) )
}


## Beta Binomial
rvalue.agrid.bb <- function( binomial.count, binomial.size, hypers,alpha.grid,
                             smooth,tst)
{
  ##############################################################
  ## computes `R' values for the conjugate beta-binomial problem
  ## via the `alpha grid' algorithm
  ## uses MLE for hypers
  ##
  ###############################################################
  ## Input
  ##
  ## binomial.count	vector of binomial counts
  ## binomial.size 	matched vector of total trials
  
  xx <- binomial.count
  nn <- binomial.size
  
  ################################################################
  ## input checks: TBD
  
  ################################################################
  ##
  ## of the conjugate beta hyper-parameters aa, bb 
  ## use values input, unless hypers="MLE"
  
  if( length(hypers)==2 & min(hypers)> 0 )
  { aa <- hypers[1]; bb <- hypers[2] }  
  if( hypers[1]=="estimate" )
  {	
    ## get MLE
    gethyps <- HyperEstimate(xx,nn,family="binomial")
    aa <- gethyps[1]
    bb <- gethyps[2]
  }
  
  # aside
  thetaPM <- ( xx+aa )/(nn+aa+bb)
  
  ################################################################
  ## alpha grid algorithm 
  ##
  nunits <- length(xx)
  ngrid <- length(alpha.grid)
  
  thetaAlpha <- qbeta( alpha.grid, shape1=aa, shape2=bb, lower.tail=FALSE )
  
  V <- matrix(NA, nunits, ngrid )
  #s1 <- aa + xx
  #s2 <- bb + (nn - xx)
  #V <- matrix(pbeta(rep(thetaAlpha,each=nunits), shape1=rep(s1,ngrid), shape2=rep(s2,ngrid), 
  #                 lower.tail=FALSE),nrow=nunits,ncol=ngrid)
  if(nunits <= ngrid) {
  ## row-wise...
      for( i in 1:nunits)
      {
         s1 <- aa + xx[i]
         s2 <- bb + (nn[i] - xx[i])
         V[i,] <- pbeta( thetaAlpha, shape1=s1, shape2=s2, lower.tail=FALSE )
      }
  }
  else {
      ## column-wise
      s1 <- aa + xx
      s2 <- bb + (nn - xx)
      for(j in 1:ngrid)  {
           V[,j] <- pbeta(thetaAlpha[j],shape1=s1,shape2=s2,lower.tail=FALSE)
      }
  }
  tmp <- rvalueGuts(dat=cbind(xx,nn), alpha.grid=alpha.grid, V=V, vfun=vfun.bb, 
                    hypers=c(aa,bb),smooth=smooth)
  rvals <- tmp$rvals
  lamfun <- tmp$lamfun
  bar <- data.frame( RValue=rvals, RV.rank=rank(rvals),
                     MLE.rank=rank(-xx/nn),
                     PM.rank=rank(-thetaPM),
                     xx=xx, nn=nn,
                     PostMean=thetaPM )
  
  ord <- order( rvals, -thetaPM )
  res <- bar[ord,]
  
  theta.quantiles <- qbeta(alpha.grid,shape1=aa,shape2=bb,lower.tail=FALSE)
  aux <- list( V=V, alpha.grid=alpha.grid,Vmarginals=tmp$lamfun,Vmarginals.smooth=tmp$smoothlamfun,
               unsorted = bar, hypers=c(aa,bb), theta.quantiles=theta.quantiles,
               prior = "conjugate", family = "binomial",smooth=smooth)
  return( list( main=res, aux=aux, rvalues=rvals) )
}


rvalue.agrid.pg <- function( poisson.count, mean.mult,hypers,alpha.grid,smooth)
{
  ##############################################################
  ## computes `R' values for the conjugate Poisson - Gamma problem
  ## via the `alpha grid' algorithm
  ##
  ###############################################################
  ## Input
  ##
  ## poisson.count	vector of binomial counts
  ## mean.mult     	matched vector of multipliers of par of interest
  
  xx <- poisson.count  ## all non-negative integers
  eta <- mean.mult  ## all positive
  
  ## model is xx[i] ~ Poisson( theta[i] * eta[i] )
  
  ################################################################
  ## input checks: TBD
  
  ################################################################
  ##
  ## use marginal max likelihood unless values are input
  
  ## of the conjugate Gamma hyper-parameters aa, bb 
  if( length(hypers)==2 & min(hypers)> 0 )
  { aa <- hypers[1]; bb <- hypers[2] }  
  if( hypers[1]=="estimate" )
  {
    ## get MLE
    gethyps <- HyperEstimate(xx,eta,family="poisson")
    aa <- gethyps[1]
    bb <- gethyps[2]
  }
  
  # aside
  thetaPM <- ( xx+aa )/(eta + bb)
  
  ################################################################
  ##
  ## alpha grid algorithm 
  ##
  nunits <- length(xx)
  ngrid <- length(alpha.grid)
  
  thetaAlpha <- qgamma( alpha.grid, shape=aa, rate=bb, lower.tail=FALSE )
  
  V <- matrix(NA, nunits, ngrid )
  ## row-wise...
  if(nunits <= ngrid) {
      for( i in 1:nunits) {
         ss <- aa + xx[i]
         rr <- bb + eta[i] 
         V[i,] <- pgamma( thetaAlpha, shape=ss, rate=rr, lower.tail=FALSE )
      }
  }
  else {
      for(j in 1:ngrid) {
          ss <- aa + xx
          rr <- bb + eta
          V[,j] <- pgamma( thetaAlpha[j], shape=ss, rate=rr, lower.tail=FALSE )
      }
  }
  
  tmp <- rvalueGuts( dat=cbind(xx,eta), alpha.grid=alpha.grid, V=V, vfun=vfun.pg, 
                     hypers=c(aa,bb),smooth=smooth)
  rvals <- tmp$rvals
  lamfun <- tmp$lamfun
  bar <- data.frame( RValue=rvals, RV.rank=rank(rvals),
                     MLE.rank=rank(-xx/eta),
                     PM.rank=rank(-thetaPM),
                     xx=xx, eta=eta,
                     PostMean=thetaPM )
  
  ord <- order( rvals, -thetaPM )
  res <- bar[ord,]
  #theta.quantiles <- qgamma(alpha.grid, shape = aa, rate = bb, lower.tail=FALSE)
  aux <- list( V=V,  alpha.grid=alpha.grid, Vmarginals=tmp$lamfun,Vmarginals.smooth=tmp$smoothlamfun,
               unsorted = bar, hypers=c(aa,bb), theta.quantiles=thetaAlpha,
               prior = "conjugate",family="poisson",smooth=smooth)
  return( list( main=res, aux=aux, rvalues=rvals ) )
}

rvalue.agrid.gg <- function(x, shapes, hypers, alpha.grid, smooth) {
  
  if( length(hypers)==2 & min(hypers) > 0 )  { 
     aa <- hypers[1] 
     bb <- hypers[2] 
  }  
  if( hypers[1]=="estimate" ) {
     ## get MLE
     gethyps <- HyperEstimate(x, shapes, family="Gamma")
     aa <- gethyps[1]
     bb <- gethyps[2]
  }
  ## Posterior Mean (Need alpha_0 + alpha > 1)
  thetaPM <- bb/((bb*x + 1)*(aa + shapes - 1))
  
  nunits <- length(x)
  ngrid <- length(alpha.grid)
  
  thetaAlpha.inv <- qgamma(alpha.grid, shape=aa, scale=bb)
  
  V <- matrix(NA, nunits, ngrid )
  ## row-wise...
  if(nunits <= ngrid) {
      for( i in 1:nunits) {
         ss <- aa + shapes[i]
         rr <- bb + x[i] 
         V[i,] <- pgamma( thetaAlpha.inv, shape=ss, rate = rr)
      }
  }
  else {
      for(j in 1:ngrid) {
          ss <- aa + shapes
          rr <- bb + x
          V[,j] <- pgamma( thetaAlpha.inv[j], shape=ss, rate = rr)
      }
  }
  tmp <- rvalueGuts(dat=cbind(x,shapes), alpha.grid=alpha.grid, V=V, vfun=vfun.gg, 
                    hypers=c(aa,bb),smooth=smooth)
  rvals <- tmp$rvals
  lamfun <- tmp$lamfun
  bar <- data.frame(RValue=rvals, RV.rank=rank(rvals),
                    MLE.rank=rank(-x/shapes), PM.rank=rank(-thetaPM),
                    xx=x, shapes=shapes, PostMean=thetaPM )
  
  ord <- order( rvals, -thetaPM )
  res <- bar[ord,]
  aux <- list(V=V, alpha.grid=alpha.grid, Vmarginals=tmp$lamfun,Vmarginals.smooth=tmp$smoothlamfun,
               unsorted = bar, hypers=c(aa,bb), theta.quantiles=1/thetaAlpha.inv,
               prior = "conjugate",family="Gamma",smooth=smooth)
  return( list( main=res, aux=aux, rvalues=rvals ) )
}




