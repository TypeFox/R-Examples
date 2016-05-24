


##
#### calculating Sigma in DPM-MVN specification
##
calVar_DPM_Normal <- function(results)
{
    J <- length(results$accept.V)
    nSample <- dim(results$class.p)[1]
    
    sigV <- rep(0,nSample)
    
    for(i in 1:nSample)
    {
        u <- length(unique(results$class.p[i,]))
        uniqC <- sort(unique(results$class.p[i,]))
        
        weight <- as.vector(table(results$class.p[i,])/J)
        mu <- rep(list(NULL), u)
        Sigma <- rep(list(NULL), u)
        
        for(j in 1:u)
        {
            mu[[j]] <- results$mu.p[i, j]
            Sigma[[j]] <- 1/results$zeta.p[i,j]
        }
        
        mu_mean <- 0
        for(j in 1:u)
        {
            mu_mean <- mu_mean + weight[j] * mu[[j]]
        }
        
        for(j in 1:u)
        {
            sigV[i] <- sigV[i] + weight[j] * ( (mu[[j]] - mu_mean)^2  + Sigma[[j]])
        }
    }
    return(sigV)
}



##
#### calculating Sigma in DPM-Normal specification
##
calVar_DPM_MVN <- function(results)
{
    J <- length(results$accept.V)
    nSample <- dim(results$class.p)[1]
    
    Sigma_V <- array(0, c(3,3,nSample))
    
    for(i in 1:nSample)
    {
        u <- length(unique(results$class.p[i,]))
        uniqC <- sort(unique(results$class.p[i,]))
        
        weight <- as.vector(table(results$class.p[i,])/J)
        mu <- rep(list(NULL), u)
        Sigma <- rep(list(NULL), u)
        
        for(j in 1:u)
        {
            mu[[j]] <- results$mu.p[i, (3*j-2):(3*j)]
            Sigma[[j]] <- results$Sigma.p[,(3*j-2):(3*j),i]
        }
        
        mu_mean <- rep(0, 3)
        for(j in 1:u)
        {
            mu_mean <- mu_mean + weight[j] * mu[[j]]
        }
        
        for(j in 1:u)
        {
            Sigma_V[,,i] <- Sigma_V[,,i] + weight[j] * ( (mu[[j]] - mu_mean) %*% t((mu[[j]] - mu_mean)) + Sigma[[j]])
        }
    }
    return(Sigma_V)
}



##
#### Potential Scale Reduction
##
calcPSR <- function(ChainMat)
{
	## number of scans
	n <- dim(ChainMat)[1]
	
	## Between and Within-chain variation
	B <- n * var(apply(ChainMat, 2, mean))
	W <- mean(apply(ChainMat, 2, var))

	## MPV = marginal posterior variance
	MPV <- ((n-1)*W + B) / n
	
	## psr = potential scale reduction
	PSR <- sqrt(MPV / W)
	return(PSR)
}




##
logLike.weibull.Uni <- function(para, y, delta, Xmat)
{
  ##
  kappa    <- exp(para[1])
  alpha <- exp(para[2])
  ##
  nP  <- length(para)
  eta <- as.vector(Xmat %*% para[3:nP])
  ##
  comp1 <- - (kappa*y^alpha) * exp(eta)
  comp2 <- log(alpha) + log(kappa) + (alpha-1) * log(y) + eta
  ##
  loglh <- sum(comp1) + sum(delta * comp2)
  ##
  return(-loglh)
}

##
logLike.weibull.SCR <- function(para, y1, y2, delta1, delta2, Xmat1=NULL, Xmat2=NULL, Xmat3=NULL, frailty=TRUE)
{
  ##
  kappa1    <- exp(para[1])
  alpha1 <- exp(para[2])
  kappa2    <- exp(para[3])
  alpha2 <- exp(para[4])  
  kappa3    <- exp(para[5])
  alpha3 <- exp(para[6])
  if(frailty == TRUE){
    theta    <- exp(para[7])
    thetaInv <- 1 / theta    
  }
  ##
  nP.0 <- ifelse(frailty, 7, 6)
  nP.1 <- ncol(Xmat1)
  nP.2 <- ncol(Xmat2)
  nP.3 <- ncol(Xmat3)
  ##
  eta.1 <- as.vector(Xmat1 %*% para[nP.0 + c(1:nP.1)])
  eta.2 <- as.vector(Xmat2 %*% para[nP.0 + nP.1 + c(1:nP.2)])
  eta.3 <- as.vector(Xmat3 %*% para[nP.0 + nP.1 + nP.2 + c(1:nP.3)])
  
  ##
  type1 <- as.numeric(delta1 == 1 & delta2 == 1)
  type2 <- as.numeric(delta1 == 0 & delta2 == 1)
  type3 <- as.numeric(delta1 == 1 & delta2 == 0)
  type4 <- as.numeric(delta1 == 0 & delta2 == 0)
  ##
  log.h1star.y1 <- log(alpha1) + log(kappa1) + (alpha1 - 1) * log(y1) + eta.1
  log.h2star.y2 <- log(alpha2) + log(kappa2) + (alpha2 - 1) * log(y2) + eta.2
  log.h3star.y2 <- log(alpha3) + log(kappa3) + (alpha3 - 1) * log(y2) + eta.3
  ##
  q.y1 <- kappa1*(y1)^alpha1 * exp(eta.1) + kappa2*(y1)^alpha2 * exp(eta.2)
  q.y2 <- kappa1*(y2)^alpha1 * exp(eta.1) + kappa2*(y2)^alpha2 * exp(eta.2)
  ##
  w.y1.y2 <- kappa3*(y2)^alpha3 * exp(eta.3) - kappa3*(y1)^alpha3 * exp(eta.3)
  ##
  if(frailty == TRUE)
  {
    logLike1 <- log.h1star.y1 + log.h3star.y2 + log(1+theta) - ((thetaInv + 2) * log(1 + (theta * (q.y1 + w.y1.y2))))
    logLike2 <- log.h2star.y2 - ((thetaInv + 1) * log(1 + (theta * q.y2)))
    logLike3 <- log.h1star.y1 - ((thetaInv + 1) * log(1 + (theta * (q.y1 + w.y1.y2))))
    logLike4 <- - thetaInv * log(1 + (theta * q.y2))
  }
  if(frailty == FALSE)
  {
    logLike1 <- log.h1star.y1 + log.h3star.y2 - (q.y1 + w.y1.y2)
    logLike2 <- log.h2star.y2 - q.y2
    logLike3 <- log.h1star.y1 - (q.y1 + w.y1.y2)
    logLike4 <- - q.y2
  }
  ##
  loglh <- sum(type1 * logLike1) + sum(type2 * logLike2) + sum(type3 * logLike3) + sum(type4 * logLike4)
  ##
  return(-loglh)
}




##
logLike.weibull.SCR.SM <- function(para, y1, y2, delta1, delta2, Xmat1=NULL, Xmat2=NULL, Xmat3=NULL, frailty=TRUE)
{
    ##
    kappa1    <- exp(para[1])
    alpha1 <- exp(para[2])
    kappa2    <- exp(para[3])
    alpha2 <- exp(para[4])
    kappa3    <- exp(para[5])
    alpha3 <- exp(para[6])
    if(frailty == TRUE){
        theta    <- exp(para[7])
        thetaInv <- 1 / theta
    }
    ##
    nP.0 <- ifelse(frailty, 7, 6)
    nP.1 <- ncol(Xmat1)
    nP.2 <- ncol(Xmat2)
    nP.3 <- ncol(Xmat3)
    ##
    eta.1 <- as.vector(Xmat1 %*% para[nP.0 + c(1:nP.1)])
    eta.2 <- as.vector(Xmat2 %*% para[nP.0 + nP.1 + c(1:nP.2)])
    eta.3 <- as.vector(Xmat3 %*% para[nP.0 + nP.1 + nP.2 + c(1:nP.3)])
    
    ##
    type1 <- as.numeric(delta1 == 1 & delta2 == 1)
    type2 <- as.numeric(delta1 == 0 & delta2 == 1)
    type3 <- as.numeric(delta1 == 1 & delta2 == 0)
    type4 <- as.numeric(delta1 == 0 & delta2 == 0)
    ##
    log.h1star.y1 <- log(alpha1) + log(kappa1) + (alpha1 - 1) * log(y1) + eta.1
    log.h2star.y2 <- log(alpha2) + log(kappa2) + (alpha2 - 1) * log(y2) + eta.2
    log.h3star.y2 <- log(alpha3) + log(kappa3) + (alpha3 - 1) * log(y2-y1) + eta.3
    ##
    q.y1 <- kappa1*(y1)^alpha1 * exp(eta.1) + kappa2*(y1)^alpha2 * exp(eta.2)
    q.y2 <- kappa1*(y2)^alpha1 * exp(eta.1) + kappa2*(y2)^alpha2 * exp(eta.2)
    ##
    w.y1.y2 <- kappa3*(y2-y1)^alpha3 * exp(eta.3)
    ##
    if(frailty == TRUE)
    {
        logLike1 <- log.h1star.y1 + log.h3star.y2 + log(1+theta) - ((thetaInv + 2) * log(1 + (theta * (q.y1 + w.y1.y2))))
        logLike2 <- log.h2star.y2 - ((thetaInv + 1) * log(1 + (theta * q.y2)))
        logLike3 <- log.h1star.y1 - ((thetaInv + 1) * log(1 + (theta * (q.y1 + w.y1.y2))))
        logLike4 <- - thetaInv * log(1 + (theta * q.y2))
    }
    if(frailty == FALSE)
    {
        logLike1 <- log.h1star.y1 + log.h3star.y2 - (q.y1 + w.y1.y2)
        logLike2 <- log.h2star.y2 - q.y2
        logLike3 <- log.h1star.y1 - (q.y1 + w.y1.y2)
        logLike4 <- - q.y2
    }
    ##
    loglh <- sum(logLike1[type1==1]) + sum(logLike2[type2==1]) + sum(logLike3[type3==1]) + sum(logLike4[type4==1])
    ##
    return(-loglh)
}

