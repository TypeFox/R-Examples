new.ctmm <- methods::setClass("ctmm", representation(info="list"), contains="list")
new.covm <- methods::setClass("covm", representation(par="numeric",isotropic="logical"), contains="matrix")
# try putting contains as unnamed first element of representation?

#######################################
# convenience wrapper for new.ctmm
ctmm <- function(tau=NULL,isotropic=FALSE,range=TRUE,circle=FALSE,CPF=FALSE,error=FALSE,...)
{
  tau.names <- c("position","velocity","acceleration")
  dim.names <- c("x","y")
  List <- list(...)
  List <- List[!sapply(List,is.null)]
  
  info <- List$info
  List$info <- NULL
  if(is.null(info)) { info=list() }
  
  List$error <- error
  List$circle <- circle
  
  # put covariance into universal format
  if(!is.null(List$sigma)) { List$sigma <- covm(List$sigma,isotropic=isotropic) }
  List$isotropic <- isotropic
  
  # label spatial elements
  # FIX THIS
  #if(!is.null(List$mu)) { List$mu <- as.numeric(List$mu) ; names(List$mu) <- dim.names }
  #if(!is.null(List$COV.mu)) { dimnames(List$COV.mu) <- list(dim.names,dim.names) }
  
  # label tau elements
  K <- length(tau)
  if(K && !CPF)
  {
    tau <- sort(tau,decreasing=TRUE)
    names(tau) <- tau.names[1:K]
    tau <- tau[tau>0]
  }
  if(!length(tau)) { tau <- NULL } # NULL, NA, integer(0), ...
  List$CPF <- CPF
  if(CPF && length(tau)) { names(tau) <- c("period","decay") }
  List$tau <- tau
  
  # label parameter covariance matrix
  cov.names <- c("area",names(tau))
  if(circle) { cov.names <- c(cov.names,"circle") }
  # List$COV can resolve to List$COV.mu apparently
  if(!is.null(List[["COV"]])) { dimnames(List$COV) <- list(cov.names,cov.names) }
  
  if(length(tau))
  {
    if(tau[1]==Inf)
    { range <- FALSE }
    else
    { range <- TRUE }
  }  
  List$range <- range

  if(!range & circle) { stop("Inconsistent model options: range=FALSE, circle=TRUE.") }
  
  # default mean function
  if(is.null(List$mean)) { List$mean <- "stationary" }
  else if(List$mean=="periodic")
  {
    P <- attr(List$mean,"par")$P
    if(is.null(P)) { attr(List$mean,"par")$P <- 24*60^2 }
    # default k via Nyquist needs data
  }
  
  # name mean parameters
  if(!is.null(List$mu))
  {
    # format simple-style stationary mean properly
    List$mu <- rbind(List$mu)
    colnames(List$mu) <- c("x","y")
    
    if(List$mean=="periodic")
    {
       k <- (nrow(List$mu)-1)/2
       K <- 1:k
       K <- c( paste("cos",K) , paste("sin",K) )
       K <- K[order(sequence(c(k,k)))]
       K <- c("const",K)
       rownames(List$mu) <- K
    }
    
    if(List$mean=="polynomial")
    {
      k <- nrow(List$mu)-1
      K <- 0:k
      rownames(List$mu) <- K  
    }

    # List$COV can resolve to List$COV.mu apparently
    if(!is.null(List[["COV.mu"]]))
    {
      if(length(dim(List$DOF.mu))==2)
      { dimnames(List$DOF.mu) <- list(dim.names,dim.names) }
      else if(length(dim(List$DOF.mu))==4)
      { dimnames(List$DOF.mu) <- list(dim.names,K,K,dim.names) }
      
      if(length(dim(List$COV.mu))==2)
      { dimnames(List$COV.mu) <- list(dim.names,dim.names) }
      else
      { dimnames(List$COV.mu) <- list(dim.names,K,K,dim.names) }
    }
  }
  
  result <- new.ctmm(List,info=info)
  
  return(result)
}
#print.ctmm <- function(x,...) { print.listof(x,...) }

ctmm.ctmm <- function(CTMM)
{
  List <- methods::getDataPart(CTMM)
  names(List) <- names(CTMM) # bug workaround
  List$info <- attr(CTMM,"info")
  CTMM <- do.call(ctmm,List)
  return(CTMM)
}

# 2D covariance matrix universal format
covm <- function(pars,isotropic=FALSE)
{
  if(is.null(pars))
  { return(NULL) }
  else if(class(pars)=="covm")
  { return(pars) }
  else if(length(pars)==1)
  {
    pars <- c(pars,0,0)
    sigma <- diag(pars[1],2)
  }
  else if(length(pars)==3)
  { sigma <- sigma.construct(pars) }
  else if(length(pars)==4)
  {
    sigma <- pars
    pars <- sigma.destruct(sigma)
  }
  
  # isotropic error check
  if(isotropic)
  {
    pars <- c(mean(diag(sigma)),0,0)
    sigma <- diag(pars[1],2)
  }
  
  name <- c("x","y")
  dimnames(sigma) <- list(name,name)
  
  names(pars) <- c("area","eccentricity","angle")
  
  new.covm(sigma,par=pars,isotropic=isotropic)
}

# construct covariance matrix from 1-3 parameters
sigma.construct <- function(pars)
{
  GM <- pars[1]
  if(length(pars)==1)
  {
    e <- 0
    theta <- 0
  }
  else
  {
    e <- pars[2]
    theta <- pars[3]
  }
  
  u <- c(cos(theta),sin(theta))
  v <- c(-sin(theta),cos(theta))
  e <- exp(e/2)
  
  sigma <- GM * ( (u%o%u)*e + (v%o%v)/e )
  
  return(sigma)
} 

# reduce covariance matrix to 1-3 parameters
sigma.destruct <- function(sigma)
{
  stuff <- eigen(sigma)
  
  e <- stuff$values
  GM <- sqrt(prod(e))
  e <- log(e[1]/e[2])
  
  if(e==0)
  { theta <- 0 }
  else
  {
    theta <- stuff$vectors[,1]
    theta <- atan(theta[2]/theta[1])
  } 
  
  return(c(GM,e,theta))
}

# blank generic function
pars <- function(...) { return(NA) }

# return the canonical parameters of a covariance matrix
pars.covm <- function(COVM)
{
  if(COVM@isotropic)
  { return(COVM@par[1]) }
  else
  { return(COVM@par) }
}

# returns the canonical parameters of a tau vector
pars.tauv <- function(tau)
{
  if(length(tau)==0)
  { return(tau) }
  else if(tau[1] < Inf)
  { return(tau) }
  else if(length(tau)==1)
  { return(NULL) }
  else
  { return(tau[-1]) }
}

############################
# coarce infinite parameters into finite parameters appropriate for numerics
###########################
ctmm.prepare <- function(data,CTMM)
{
  K <- length(CTMM$tau)  # dimension of hidden state per spatial dimension
  
  range <- TRUE
  
  if(K>0)
  { 
    # numerical limit for rangeless processes
    if(CTMM$tau[1]==Inf)
    {
      range <- FALSE
      CTMM$tau[1] <- log(2^(52/2))*(last(data$t)-data$t[1])
      
      # diffusion -> variance
      if(!is.null(CTMM$sigma))
      { 
        T <- (CTMM$tau[1]-if(K>1){CTMM$tau[2]}else{0})
        CTMM$sigma <- CTMM$sigma*T
        CTMM$sigma@par[1] <- CTMM$sigma@par[1]*T
      }
    }
    
    # continuity reduction
    CTMM$tau = CTMM$tau[CTMM$tau!=0]
    K <- length(CTMM$tau)
  }
  # I am this lazy
  if(K==0) { K <- 1 ; CTMM$tau <- 0 }
  
  CTMM$range <- range
  
  # evaluate mean function for this data set if no vector is provided
  if(is.null(CTMM$mean.vec))
  { CTMM$mean.vec <- call.mean(CTMM,data$t) }
  
  return(CTMM)
}

ctmm.repair <- function(CTMM)
{
  if(!CTMM$range)
  {
    K <- length(CTMM$tau)
    
    # variance -> diffusion
    T <- (CTMM$tau[1]-if(K>1){CTMM$tau[2]}else{0})
    CTMM$sigma <- CTMM$sigma/T
    CTMM$sigma@par[1] <- CTMM$sigma@par[1]/T
    
    CTMM$tau[1] <- Inf
    
    # delete garbate estimates
    CTMM$mu <- NULL
    CTMM$COV.mu <- NULL
    CTMM$DOF.mu <- NULL
  }
  
  # erase evaluated mean vector from ctmm.prepare
  CTMM$mean.vec <- NULL
  
  return(CTMM)
}

## prepare error array
error.prepare <- function(DATA,CTMM)
{
  n <- length(DATA$t)
  
  # model the error
  if(CTMM$error>0)
  {
    # is the data supplied with error estimates
    error <- DATA$error
    # if not, then use the modeled error
    if(is.null(error)) { error <- rep(as.numeric(CTMM$error),n) }
  }
  else
  { error <- rep(0,n) }
  
  return(error)
}

# degree of continuity in the model
continuity <- function(CTMM)
{
  K <- sum(CTMM$tau > 0)
  return(K)
}

###############################
# Propagator/Green's function and Two-time correlation from Langevin equation for Kalman filter and simulations
langevin <- function(dt,CTMM)
{
  tau <- CTMM$tau
  CPF <- CTMM$CPF
  sigma <- methods::getDataPart(CTMM$sigma)
  K <- continuity(CTMM)
  
  if(K==0)
  {
    Green <- array(0,c(1,1))
    Sigma <- array(1,c(1,1))
  }
  else if(K==1)
  {
    c0 <- exp(-dt/tau)
    Green <- array(c0,c(1,1))
    Sigma <- array(1-c0^2,c(1,1))
  }
  else if(K==2)
  {
    f <- 1/tau
    
    if(CPF)
    {
      nu <- 2*pi*f[1]
      f <- f[2]
      Omega2 <- f^2 + nu^2
    }
    else
    { Omega2 <- 1/prod(tau) }

    T <- 2*mean(f)/Omega2
    
    if(dt==Inf) # make this numerically relative in future
    {
      Green <- rbind( c(0,0) , c(0,0) )
      Sigma <- rbind( c(1,0) , c(0,Omega2) )
    }
    else
    {
      if(CPF)
      {
        Exp <- exp(-f*dt)
        SinE <- sin(nu*dt)*Exp
        CosE <- cos(nu*dt)*Exp
        
        c0 <- CosE + (f/nu)*SinE
        c1 <- -(Omega2/nu)*SinE
        c2 <- -Omega2*(CosE - (f/nu)*SinE)
      }
      else
      {
        Exp <- exp(-dt/tau)/diff(tau)
        c0 <- diff(Exp*tau)
        c1 <- -diff(Exp)
        c2 <- diff(Exp/tau)
      }
      
      Green <- rbind( c(c0,-c1/Omega2) , c(c1,-c2/Omega2) )
      Sigma <- -T*c1^2  #off-diagonal term
      Sigma <- rbind( c(1,0)-c(c0^2+c1^2/Omega2,Sigma) , c(0,Omega2)-c(Sigma,c1^2+c2^2/Omega2) )
    }
  }
  
  return(list(Green=Green, Sigma=sigma*Sigma))
}

#############################################################
# Internal Kalman filter/smoother for multiple derivatives, dimensions, trends
# Kalman filter/smoother for matrix P operator and multiple mean functions
# this is a lower level function
# more for generalizability/readability than speed at this point
kalman <- function(z,u,dt,CTMM,error=rep(0,nrow(z)),smooth=FALSE,sample=FALSE)
{
  n <- nrow(z)
  DATA <- 1:ncol(z)
  
  # glob data and mean functions together for Kalman filter
  z <- cbind(z,u)
  VEC <- ncol(z)
  
  # indices of mean functions
  MEAN <- (last(DATA)+1):VEC

  tau <- CTMM$tau
  K <- length(tau)  # dimension of hidden state per spatial dimension

  # observed dimensions
  OBS <- 1
  Id <- diag(OBS)
  
  # observable state projection operator (will eventially upgrade to use velocity telemetry)
  P <- array(0,c(K,OBS))
  P[1:OBS,] <- 1
  
  # forecast estimates for zero-mean, unit-variance filters
  zFor <- array(0,c(n,K,VEC))
  sFor <- array(0,c(n,K,K))
  
  # forcast residuals
  zRes <- array(0,c(n,OBS,VEC))
  sRes <- array(0,c(n,OBS,OBS))
  
  # concurrent/filtered estimates
  zCon <- array(0,c(n,K,VEC))
  sCon <- array(0,c(n,K,K))
  
  # initial state info
  Langevin <- langevin(dt=dt[1],CTMM=CTMM)
  Green <- Langevin$Green
  Sigma <- Langevin$Sigma
  sFor[1,,] <- Sigma
  
  for(i in 1:n)
  {
    # residual covariance
    sForP <- sFor[i,,] %*% P # why do I need this?
    sRes[i,,] <- ((t(P) %*% sForP) + error[i]*Id)
    
    # forcast residuals
    zRes[i,,] <- z[i,] - (t(P) %*% zFor[i,,])
    
    if(all(abs(sRes[i,,])<Inf)){ Gain <- sForP %*% solve(sRes[i,,]) }
    else { Gain <- sForP %*% (0*Id) } # solve() doesn't like this case...
    
    # concurrent estimates
    zCon[i,,] <- zFor[i,,] + (Gain %*% zRes[i,,])
    sCon[i,,] <- (sFor[i,,] - (Gain %*% t(sForP)))
    
    # update forcast estimates for next iteration
    if(i<n)
    {
      # does the time lag change values? Then update the propagators.
      if(dt[i+1] != dt[i])
      {
        Langevin <- langevin(dt=dt[i+1],CTMM=CTMM)
        Green <- Langevin$Green
        Sigma <- Langevin$Sigma
      }
      #update forcast estimates now
      zFor[i+1,,] <- Green %*% zCon[i,,]
      sFor[i+1,,] <- ((Green %*% sCon[i,,] %*% t(Green)) + Sigma)
    }
  }
  
  # general quadratic form, tracing over times
  # hard-code weights for location observation only
  M <- Adj(zRes[,1,]) %*% (zRes[,1,]/sRes[,1,1])

  # estimate mean parameter
  W <- as.matrix(M[MEAN,MEAN])
  mu <- solve(W) %*% M[MEAN,DATA]

  # returned profiled mean
  if(!smooth && !sample)
  {
    # hard coded for position observations
    sigma <- (M[DATA,DATA] - (Adj(mu) %*% W %*% mu))/n

    # log det autocorrelation matrix == trace log autocorrelation matrix
    logdet <- sum(log(sRes))
    
    return(list(mu=mu,W=W,sigma=sigma,logdet=logdet))
  }
  # delete residuals
  rm(zRes,sRes)
  
  #####################
  # KALMAN SMOOTHER
  #####################
  # Finish detrending the effect of a stationary mean
  MU <- zFor[,,MEAN]
  dim(MU) <- c(n*K,length(MEAN))
  MU <- MU %*% mu
  dim(MU) <- c(n,K,length(DATA))
  zFor[,,DATA] <- zFor[,,DATA,drop=FALSE] - MU
  # there has to be a better way to do this?
  MU <- zCon[,,MEAN]
  dim(MU) <- c(n*K,length(MEAN))
  MU <- MU %*% mu
  dim(MU) <- c(n,K,length(DATA))
  zCon[,,DATA] <- zCon[,,DATA,drop=FALSE] - MU
  # why does R drop dimensions so randomly?

  # delete u(t)
  zCon <- zCon[,,DATA,drop=FALSE]
  zFor <- zFor[,,DATA,drop=FALSE]
  # drop=FALSE must be here for BM/OU and I don't fully understand why
  
  # upgrade concurrent estimates to Kriged estimates  
  Green <- langevin(dt=dt[n],CTMM=CTMM)$Green
  for(i in (n-1):1)
  {
    L <- sCon[i,,] %*% t(Green) %*% solve(sFor[i+1,,])
    
    # overwrite concurrent estimate with smoothed estimate
    zCon[i,,] <- zCon[i,,] + L %*% (zCon[i+1,,]-zFor[i+1,,])
    sCon[i,,] <- sCon[i,,] + L %*% (sCon[i+1,,]-sFor[i+1,,]) %*% t(L)
    
    # next time's propagator if necessary
    if(dt[i] != dt[i+1])
    { Green <- langevin(dt=dt[i],CTMM=CTMM)$Green }
    
    #################
    # RANDOM SAMPLER
    #################
    if(sample)
    {
      zCon[i,,] <- sapply(DATA,function(d){MASS::mvrnorm(mu=zCon[i,,d],Sigma=sCon[i,,])})
      sCon[i,,] <- array(0,c(K,K))
    }
  }
  
  # restore stationary mean to locations only
  zCon[,1,] <- zCon[,1,] + (cbind(u) %*% mu)
  
  zname <- c("position")
  if(K>1) { zname <- c(zname,"velocity") }
  dimnames(zCon) <- list(NULL,zname,c("x","y")[DATA])
  dimnames(sCon) <- list(NULL,zname,zname)
  
  # return smoothed states
  # this object is temporary
  state <- list(CTMM=CTMM,Z=zCon,S=sCon)
  class(state) <- "state"
  
  return(state)
}

####################################
# log likelihood function
####################################
ctmm.loglike <- function(data,CTMM=ctmm(),verbose=FALSE)
{
  n <- length(data$t)
  
  # prepare model for numerics
  CTMM <- ctmm.prepare(data,CTMM)
  
  range <- CTMM$range
  isotropic <- CTMM$isotropic

  sigma <- CTMM$sigma
  if(!is.null(sigma))
  {
    area <- sigma@par[1]
    ecc <- sigma@par[2]
    theta <- sigma@par[3]
  }
  else
  {
    area <- NA
    ecc <- 0
    theta <- 0
  }
  
  circle <- CTMM$circle
  if(circle) { circle <- 2*pi/circle }
  if(abs(circle) == Inf) { circle <- FALSE }

  n <- length(data$t)

  t <- data$t
  # time lags
  dt <- c(Inf,diff(t))
  
  # data z and mean vector u
  z <- cbind(data$x,data$y)
  u <- CTMM$mean.vec
  M <- ncol(u) # number of linear parameters per spatial dimension

  error <- error.prepare(data,CTMM)

  # do we need to orient the data along the major an minor axes of sigma
  ROTATE <- !isotropic && (CTMM$error || circle)
  if(ROTATE) { z <- z %*% t(rotate(-theta)) }
  
  if(circle) # ONE KALMAN FILTER WITH COMPLEX SIGNAL
  {
    # proportional standardization from ellipse to circle
    if(ecc)
    {
      z[,1] <- z[,1] * exp(-ecc/4)
      z[,2] <- z[,2] * exp(+ecc/4)
    }
    z <- cbind(z[,1] + 1i*z[,2])
    
    # corotating frame
    R <- exp(-1i*circle*(t-t[1]))
    z <- R * z
    u <- R * u
    
    if(!CTMM$error)
    {
      CTMM$sigma <- 1
      KALMAN <- kalman(z,u,dt=dt,CTMM=CTMM,error=error)
      
      ML.area <- as.numeric(Re(KALMAN$sigma))/2
      # profile variance if unspecified
      if(is.na(area))
      { 
        area <- ML.area
        sigma <- covm(c(area,ecc,theta),isotropic=isotropic)
      }
      
      COV.mu <- area * solve(KALMAN$W)
      DOF.mu <- KALMAN$W
      
      loglike <- -KALMAN$logdet -(n)*log(2*pi*area) - (n)*(ML.area/area)
    }
    else
    {
      CTMM$sigma <- area
      KALMAN <- kalman(z,u,dt=dt,CTMM=CTMM,error=error)
      
      R.sigma <- KALMAN$sigma
      
      COV.mu <- solve(KALMAN$W)
      DOF.mu <- area * KALMAN$W
      
      loglike <- -KALMAN$logdet - (n)*log(2*pi*1) - (n/2)*R.sigma
    }

    loglike <- Re(loglike)
    # real array formatting
    mu <- KALMAN$mu
    mu <- cbind(Re(mu),Im(mu))
    # complex correlations are x-y correlations
    COV.mu <- array(c(Re(COV.mu),-Im(COV.mu),Im(COV.mu),Re(COV.mu)),c(M,M,2,2))
    DOF.mu <- array(c(Re(DOF.mu),-Im(DOF.mu),Im(DOF.mu),Re(DOF.mu)),c(M,M,2,2))
    
    # de-standardization
    R <- exp(c(+1,-1)*ecc/4)
    mu <- t(R * t(mu))
    COV.mu <- array(COV.mu,c(M^2*2,2))
    COV.mu <- t(R * t(COV.mu))
    COV.mu <- array(COV.mu,c(M,M,2,2))
    COV.mu <- aperm(COV.mu,c(1,2,4,3))
    COV.mu <- array(COV.mu,c(M^2*2,2))
    COV.mu <- t(R * t(COV.mu))
    COV.mu <- array(COV.mu,c(M,M,2,2))
  }
  else if(!CTMM$error) # ONE KALMAN FILTER WITH NO ERROR
  {
    CTMM$sigma <- 1
    KALMAN <- kalman(z,u,dt=dt,CTMM=CTMM,error=error)
    
    ML.sigma <- KALMAN$sigma
    # profile covariance if sigma unspecified
    if(is.null(sigma)) { sigma <- covm(ML.sigma,isotropic=isotropic) }
    
    mu <- KALMAN$mu
    COV.mu <- solve(KALMAN$W) %o% sigma
    DOF.mu <- KALMAN$W %o% diag(2)
    
    loglike <- -KALMAN$logdet -(n/2)*log(det(2*pi*sigma)) - (n/2)*sum(diag(ML.sigma %*% solve(sigma)))
  }
  else if(isotropic) # ONE KALMAN FILTER WITH ERROR
  {
    CTMM$sigma <- area
    KALMAN <- kalman(z,u,dt=dt,CTMM=CTMM,error=error)
    
    mu <- KALMAN$mu
    COV.mu <- solve(KALMAN$W) %o% diag(2)
    DOF.mu <- (area * KALMAN$W) %o% diag(2)
    
    # includes error and standardization
    R.sigma <- KALMAN$sigma
    R.sigma <- mean(diag(R.sigma))
    
    loglike <- -KALMAN$logdet -(n)*log(2*pi*1) -(n)*R.sigma
  }
  else # TWO KALMAN FILTERS WITH ERROR
  {
    # eigen variances
    SIGMA <- area * exp(c(+1,-1)*ecc/2)
    
    # major axis likelihood
    CTMM$sigma <- SIGMA[1]
    KALMAN1 <- kalman(cbind(z[,1]),u,dt=dt,CTMM=CTMM,error=error)
    
    # minor axis likelihood
    CTMM$sigma <- SIGMA[2]
    KALMAN2 <- kalman(cbind(z[,2]),u,dt=dt,CTMM=CTMM,error=error)
    
    logdet <- KALMAN1$logdet + KALMAN2$logdet
    R.sigma <- KALMAN1$sigma + KALMAN2$sigma # standardized residual variances
    
    mu <- cbind(KALMAN1$mu,KALMAN2$mu)
    COV.mu <- array(c(solve(KALMAN1$W),diag(0,M),diag(0,M),solve(KALMAN2$W)),c(M,M,2,2)) # -1/Hessian
    DOF.mu <- array(c(SIGMA[1]*KALMAN1$W,diag(0,M),diag(0,M),SIGMA[2]*KALMAN2$W),c(M,M,2,2))

    loglike <- -(1/2)*logdet - (n)*log(2*pi*1) - (n/2)*R.sigma
  }
  
  # restructure indices from m,n,x,y to x,m,n,y
  COV.mu <- aperm( COV.mu , c(3,1,2,4))
  DOF.mu <- aperm( DOF.mu , c(3,1,2,4))
  
  # transform results back
  if(ROTATE)
  {
    R <- rotate(+theta)
    
    mu <- mu %*% t(R)
    
    COV.mu <- array(COV.mu,c(2,M^2*2))
    DOF.mu <- array(DOF.mu,c(2,M^2*2))
    
    COV.mu <- R %*% COV.mu
    DOF.mu <- R %*% DOF.mu
    
    COV.mu <- array(COV.mu,c(2*M^2,2))
    DOF.mu <- array(DOF.mu,c(2*M^2,2))

    COV.mu <- COV.mu %*% t(R)
    DOF.mu <- DOF.mu %*% t(R)
    
    COV.mu <- array(COV.mu,c(2,M,M,2))
    DOF.mu <- array(DOF.mu,c(2,M,M,2))
  }
  
  # should I drop the indices in COV.mu and DOF.mu if possible ?
  COV.mu <- drop(COV.mu)
  DOF.mu <- drop(DOF.mu)
  
  # isotropic reduction if possible
  if(length(dim(DOF.mu))==2)
  { if(DOF.mu[1,1]==DOF.mu[2,2] && DOF.mu[1,2]==0) { DOF.mu <- mean(diag(DOF.mu)) } }
    
  if(verbose)
  {
    # assign variables
    CTMM$sigma <- sigma
    CTMM <- ctmm.repair(CTMM)
    
    if(range)
    {
      CTMM$mu <- mu
      CTMM$COV.mu <- COV.mu
      CTMM$DOF.mu <- DOF.mu
    }
    
    CTMM$loglike <- loglike
    attr(CTMM,"info") <- attr(data,"info")

    CTMM <- ctmm.ctmm(CTMM)
    return(CTMM)
  }
  else  { return(loglike) }
}


###########################################################
# FIT MODEL WITH LIKELIHOOD FUNCTION (convenience wrapper to optim)
ctmm.fit <- function(data,CTMM=ctmm(),debias=TRUE,control=list(maxit=.Machine$integer.max),...)
{
  # basic info
  n <- length(data$t)
  tau <- CTMM$tau
  K <- length(tau)
  CPF <- CTMM$CPF
  circle <- CTMM$circle
  sigma <- CTMM$sigma
  if(!is.null(sigma)) { sigma <- sigma@par }
  CTMM$mu <- NULL # can always profile mu analytically
  isotropic <- CTMM$isotropic
  error <- CTMM$error
  range <- CTMM$range

  # erase previous fitting info if present
  CTMM$COV <- NULL
  CTMM$COV.mu <- NULL
  CTMM$DOF.mu <- NULL
    
  # evaluate mean function for this data set if no vector is provided
  if(is.null(CTMM$mean.vec))
  { CTMM$mean.vec <- call.mean(CTMM,data$t) }
  
  # CAVEATS
  if(error && circle && !isotropic) { warning("error==TRUE & circle==TRUE & isotropic==FALSE distorts error circle for speed.") }
  if(circle && !range) { stop("circle==TRUE & range==FALSE are incompatible.") }
  
  # parameter indices for non-profiled parameters that we have to numerically optimize
  # PARAMETERS THAT WE MAY OR MAY NOT DIFFERENTIATE WRT
  if(error) # profile 1-3 parameters
  {
    SIGMA <- 1:(if(isotropic){1}else{3})
    SIGMAV <- 1:(if(isotropic){1}else{3})
  }
  else if(circle) # profile 0-2 parameters
  {
    SIGMA <- (if(isotropic){NULL}else{1:2})
    SIGMAV <- (if(isotropic){NULL}else{2:3})
  }
  else
  { SIGMA <- NULL } 

  # PARAMETERS THAT WE WILL DIFFERENTIATE WRT
  NAMES <- NULL
  TAU <- NULL
  CIRCLE <- NULL
  calPARS <- function()
  {
    NAMES <<- "sigma"
    
    if(K+range==1) # IID, BM
    { TAU <<- NULL }
    else
    {
      NAMES <<- c(NAMES,paste("tau",names(tau)[(1+(1-range)):K]))
      TAU <<- length(SIGMA) + 1:(K-(1-range))
    }
    
    if(circle)
    {
      NAMES <<- c(NAMES,"circle")
      CIRCLE <<- length(SIGMA) + length(TAU) + 1 
    }
  }
  calPARS()

  # PARAMETERS THAT WE WILL NOT DIFFERENTIATE WRT
  # If error is an error estimate rather than TRUE, and if there is no error annotated, then we will fit error
  if(error>0 && !any(data$error>0))
  { ERROR <- length(SIGMA) + length(TAU) + length(CIRCLE) + 1 }
  else
  { ERROR <- NULL }
  
  # numerically fit parameters
  PARS <- length(SIGMA) + length(TAU) + length(CIRCLE) + length(ERROR)
  # degrees of freedom, including the mean, variance/covariance, tau, and error model
  k.mean <- ncol(CTMM$mean.vec)
  k <- 2*(k.mean - (if(range){0}else{1})) + (if(isotropic){1}else{3}) + length(TAU) + length(CIRCLE) + length(ERROR)
  
  # OPTIMIZATION GUESS (pars)
  # also construct reasonable parscale
  pars <- NULL
  parscale <- NULL
  calpars <- function()
  {
    pars <<- NULL
    parscale <<- NULL
    
    if(length(SIGMA)>0)
    {
      # need some initial guess...
      if(is.null(sigma)) { sigma <<- covm(stats::cov(cbind(data$x,data$y)),isotropic=isotropic)@par }
      pars <<- sigma[SIGMAV]
      parscale <<- c(sigma[1],1,pi/4)[SIGMAV]
      
      # can we profile the variance, if so delete the guess
      if(!error) { sigma[1] <<- NA }
    }
    
    if(length(TAU)>0)
    {
      pars <<- c(pars,pars.tauv(tau))
      parscale <<- c(parscale,pars[TAU])
    }

    # use 1/T as circulation parameter
    if(length(CIRCLE)>0)
    {
      pars <<- c(pars,1/circle)
      parscale <<- c(parscale,1/abs(circle))
    }
    
    if(length(ERROR)>0)
    {
      #assuming GPS for now
      pars <<- c(pars,error)
      parscale <<- c(parscale,error) 
    }
  }
  calpars()
  
  # OPTIMIZATION FUNCTION (fn)
  # optional argument lengths: TAU, TAU+1, TAU+SIGMA
  fn <- function(p)
  {
    # Fix sigma if provided, up to degree provided
    if(length(SIGMA)==0)
    { sigma <- NULL }
    else
    {
      # write over inherited sigma with any specified parameters
      sigma[SIGMAV] <- p[SIGMA]
      
      # enforce positivity
      sigma[-3] <- abs(sigma[-3])
      
      # for some reason optim jumps to crazy big parameters
      if(!is.na(sigma[1]) && sigma[1]*exp(abs(sigma[2])/2)==Inf) { return(Inf) }
    }
    
    # fix tau from par
    if(length(TAU)==0)
    { 
      if(range) { tau <- NULL }
      else { tau <- Inf }
    }
    else
    {
      tau <- p[TAU]
      
      # enforce positivity
      tau <- abs(tau)
      
      if(!CPF) { tau <- sort(tau,decreasing=TRUE) }
      
      if(!range) { tau <- c(Inf,tau) }
    }
    
    # fix circulation from par
    if(length(CIRCLE)>0) { circle <- 1/p[CIRCLE] }

    # fix error from par
    if(length(ERROR)>0)
    {
      error <- p[ERROR]
      
      # enforce positivity
      error <- abs(error)
    }

    # overwrite modified parameters inside this function's environment only
    CTMM$tau <- tau
    CTMM$circle <- circle
    CTMM$sigma <- covm(sigma,isotropic=isotropic)
    CTMM$error <- error
    
    return(-ctmm.loglike(data,CTMM))
  }
  
  # NOW OPTIMIZE
  if(PARS==0) # EXACT
  {
    # Bi-variate Gaussian || Brownian motion with zero error
    CTMM <- ctmm.loglike(data,CTMM=CTMM,verbose=TRUE)
    
    if(K==0){ DOF <- n ; CTMM$tau <- NULL }
    else { DOF <- n-1 }
    
    GM <- CTMM$sigma@par[1]
    CTMM$COV <- rbind(if(isotropic) { GM^2/DOF } else { GM^2/DOF/2 })
  }
  else # all further cases require optimization
  {
    if(PARS==1) # Brent is the best choice here
    {
      RESULT <- NULL
      # direct attempt that can caputure zero boundary
      ATTEMPT <- stats::optim(par=pars,fn=fn,method="Brent",lower=0,upper=10*pars)
      RESULT <- rbind(RESULT,c(ATTEMPT$par,ATTEMPT$value))
      # log scale backup that can't capture zero boundary
      ATTEMPT <- stats::nlm(function(p){f=fn(pars*exp(p))},p=0,stepmax=log(10),iterlim=control$maxit)
      RESULT <- rbind(RESULT,c(pars*exp(ATTEMPT$estimate),ATTEMPT$minimum))
      # choose the better estimate
      MIN <- sort(RESULT[,2],index.return=TRUE)$ix[1]
      pars <- RESULT[MIN,1]
    }
    else # Nelder-Mead is generally the safest and is default
    {
      control$parscale <- parscale
      pars <- stats::optim(par=pars,fn=fn,control=control,...)$par
    }
    
    # write best estimates over initial guess
    tau <- abs(pars[TAU])
    if(!range){ tau <- c(Inf,tau) }
    
    # save circulation if numerically optimized
    if(circle)
    {
      if(pars[CIRCLE])
      { circle <- 1/pars[CIRCLE] }
      else
      { circle <- FALSE }
      # In case ML circulation is zero, deactivate it in the model
    }
    
    # save sigma if numerically optimized
    if(length(SIGMA))
    {
      sigma[SIGMAV] <- pars[SIGMA]
      # enforce positivity
      sigma[-3] <- abs(sigma[-3])
    }
    else
    { sigma <- NULL }
    
    # save error magnitude if modeled
    if(length(ERROR)>0)
    { error <- abs(pars[ERROR]) }
    
    CTMM$tau <- tau
    CTMM$circle <- circle
    CTMM$sigma <- covm(sigma,isotropic=isotropic)
    CTMM$error <- error
    
    # verbose ML information
    CTMM <- ctmm.loglike(data,CTMM,verbose=TRUE)
    
    # calculate area-tau uncertainty
    sigma <- CTMM$sigma@par
    SIGMA <- 1
    SIGMAV <- 1
    calPARS()
    ERROR <- NULL
    calpars()
    
    hess <- numDeriv::hessian(fn,pars)
    grad <- numDeriv::grad(fn,pars,side=sign(pars))
    # robust covariance calculation
    CTMM$COV <- cov.loglike(hess,grad)
    # convert from circulation frequency to circulation period
    if(circle)
    { 
      g <- -circle^2
      CTMM$COV[CIRCLE,] <- g*CTMM$COV[CIRCLE,]
      CTMM$COV[,CIRCLE] <- g*CTMM$COV[,CIRCLE]
    }
    dimnames(CTMM$COV) <- list(NAMES,NAMES)
  }
  
  CTMM$mean.vec <- NULL
  
  CTMM$AICc <- (2*k-2*CTMM$loglike) + 2*k*(k+1)/(n-k-1)

  # fix for range==FALSE
  if(debias)
  {
    A <- n/(n-k.mean)
    
    sigma <- methods::getDataPart(CTMM$sigma)
    sigma <- A * sigma
    CTMM$sigma <- covm(sigma)
    
    CTMM$COV.mu <- A^2 * CTMM$COV.mu
    CTMM$COV[1,] <- A * CTMM$COV[1,]
    CTMM$COV[,1] <- A * CTMM$COV[,1]
  }
  
  return(CTMM)
}


####################################
# Newton-Raphson iterate to a ctmm model
# to test if optim worked
newton.ctmm <- function(data,CTMM)
{
  tau <- CTMM$tau
  
  # wrapper function to differentiate
  fn <- function(par)
  { 
    # will need to update this for telemetry error
    return(-ctmm.loglike(data,CTMM=ctmm(tau=par)))
  }
  
  D <- numDeriv::grad(fn,tau)
  H <- numDeriv::hessian(fn,tau)
  
  tau <- tau - (H %*% D)

  return(ctmm(tau=tau,info=attr(data,"info")))
}


###################################################
# Calculate good CIs for other functions
confint.ctmm <- function(model,alpha=0.05)
{
  # z-values for low, ML, high estimates
  z <- stats::qnorm(1-alpha/2)
  
  tau <- model$tau
  tau <- tau[tau<Inf]
  K <- length(tau)
  
  COV <- model$COV
  
  par <- NULL
  NAME <- NULL
  
  # timescale uncertainty: can hit 0 and Inf
  if(K>0)
  {
    for(k in 1:K)
    {
      NAME <- c(NAME,paste("tau",names(tau[k])))
      par <- rbind(par,confint.tau(tau[k],COV[k+1,k+1],z))
    }
  }
  
  # circulation period
  circle <- model$circle
  if(circle)
  {
    NAME <- c(NAME,"circle")
    par <- rbind(par,confint.tau(circle,COV["circle","circle"],z,min=-Inf))
  }  
    
  # standard area uncertainty: chi-square
  GM.sigma <- sqrt(det(model$sigma))
  COV <- COV[1,1] 
  # effective DOF derived from ML curvature
  DOF <- 2*GM.sigma^2/COV
  
  NAME <- c("area",NAME)
  par <- rbind(chisq.ci(GM.sigma,COV,alpha),par)
  
  rownames(par) <- NAME
  
  return(par)
}

# timescale uncertainty: can hit 0 and Inf
confint.tau <- function(tau,COV,z,min=0,max=Inf)
{
  # tau normal for lower CI
  CI <- tau + c(1,-1)*z*sqrt(COV)
  
  # lower CI of f==1/tau normal for upper CI of tau
  CI <- c(CI, (1/tau + c(1,-1)*z*sqrt(COV/tau^4))^-1)

  # take most conservative estimates
  CI <- range(CI)
  
  # enforce constraints
  CI <- c(max(CI[1],min),min(CI[2],max))
  
  CI <- c(CI[1],tau,CI[2])
  return(CI)
}

summary.ctmm <- function(object,level=0.95,level.UD=0.95,...)
{
  CLASS <- class(object)
  if(CLASS=="ctmm")
  { return(summary.ctmm.single(object,level=level,level.UD=level.UD)) }
  else if(CLASS=="list")
  { return(summary.ctmm.list(object,level=level,level.UD=level.UD)) }
}
  

######################################################
# Summarize results
summary.ctmm.single <- function(object, level=0.95, level.UD=0.95, ...)
{
  alpha <- 1-level
  alpha.UD <- 1-level
  
  # z-values for low, ML, high estimates
  z <- stats::qnorm(1-alpha/2)*c(-1,0,1)
  
  CPF <- object$CPF
  circle <- object$circle
  tau <- object$tau
  if(length(tau)>0 && tau[1]==Inf) { range <- FALSE } else { range <- TRUE }
  tau <- tau[tau<Inf]
  K <- length(tau)
  
  AM.sigma <- mean(diag(object$sigma))
  GM.sigma <- sqrt(det(object$sigma))
  ecc <- object$sigma@par[2]
  
  COV <- object$COV
  P <- nrow(COV)
  
  # where to store unit information
  name <- rep("",K+1)
  scale <- rep(1,K+1)
  
  par <- confint.ctmm(object,alpha=alpha)
  
  # standard area to home-range area
  par[1,] <- -2*log(alpha.UD)*pi*par[1,]
  
  # pretty area units
  unit.list <- unit(par[1,2],"area")
  name[1] <- unit.list$name
  scale[1] <- unit.list$scale
  
  # pretty time units
  P <- nrow(par)
  if(P>1)
  {
    for(i in 2:P)
    {
      unit.list <- unit(par[i,2],"time")
      name[i] <- unit.list$name
      scale[i] <- unit.list$scale
    }
  }
  
  # can we estimate speed?
  if(K>1 || (!range && K>0))
  {
    # RMS velocity
    if(CPF)
    {
      Omega2 <- sum(c(2*pi,1)^2/tau^2)
      grad <- -2*c(2*pi,1)^2/tau^3
    }
    else
    {
      Omega2 <- 1/prod(tau)
      grad <- -Omega2/tau
    }

    # contribution from circulation
    omega2 <- 0
    if(circle)
    {
      omega2 <- (2*pi/circle)^2
      grad <- c(grad,-2*omega2/circle)
    }
    
    # contribution from sigma
    # GM.sigma <- cosh(ecc)*AM.sigma
    ms <- AM.sigma*(Omega2+omega2)
    grad <- c(cosh(ecc)*(Omega2+omega2), AM.sigma*grad)
    var.ms <- (grad) %*% COV %*% (grad)
    # include mean
    MSPEED <- mspeed(object)
    ms <- ms + MSPEED$MS
    var.ms <- var.ms + MSPEED$VAR
    # root mean square velocity
    rms <- sqrt(ms)
    var.rms <- var.ms * (1/2/rms)^2
    
    # pretty units
    unit.list <- unit(rms,"speed")
    name <- c(name,unit.list$name)
    scale <- c(scale,unit.list$scale)
    
    log.rms <- log(rms)
    var.log.rms <- var.rms * (1/rms)^2
    
    log.rms <- log.rms + z*sqrt(var.log.rms)
    rms <- exp(log.rms)
    
    par <- rbind(par,rms)
    rownames(par)[P+1] <- "speed"
  }
  
  # Fix unit choice
  par <- par/scale
  
  # affix units
  rownames(par) <- paste(rownames(par)," (",name,")",sep="")
  
  colnames(par) <- c("low","ML","high")

  if(!range) { par <- par[-1,] } # delete off "area" (really diffusion)
  
  return(par)
}
#methods::setMethod("summary",signature(object="ctmm"), function(object,...) summary.ctmm(object,...))


sort.ctmm <- function(x, decreasing=FALSE, IC="AICc", ...)
{
  ICS <- sapply(x,function(m){m[[IC]]})
  IND <- sort(ICS,method="quick",index.return=TRUE,decreasing=decreasing)$ix
  x <- x[IND]
}


summary.ctmm.list <- function(object, IC="AICc", ...)
{
  object <- sort.ctmm(object,IC=IC)
  ICS <- sapply(object,function(m){m[[IC]]})
  ICS <- ICS - ICS[[1]]
  ICS <- array(ICS,c(length(ICS),1))
  rownames(ICS) <- names(object)
  colnames(ICS) <- paste("d",IC,sep="")
  return(ICS)
}

###################################################
# fit a bunch of models and select the best one
ctmm.select <- function(data,CTMM,verbose=FALSE,IC="AICc",...)
{
  MODELS <- list()

  # listify model guesses
  if(class(CTMM)=="ctmm") { CTMM <- list(CTMM) }
  
  for(i in 1:length(CTMM))
  {  
    if(CTMM[[i]]$CPF)
    { MODELS <- c(MODELS,ctmm.select.cpf(data,CTMM[[i]])) }
    else
    { MODELS <- c(MODELS,ctmm.select.ouf(data,CTMM[[i]])) }
  }
  
  # fit models
  MODELS <- lapply(MODELS,function(M) { ctmm.fit(data,M,...) })
  
  # sort models by AICc
  MODELS <- sort.ctmm(MODELS,IC=IC)
  
  # return everything
  if(verbose) { return(MODELS) }
  else { return(MODELS[[1]]) }
}


ctmm.select.cpf <- function(data,CTMM)
{
  MODELS <- list()
  
  MOD <- CTMM
  MODELS <- c(MODELS,ctmm.select.isotropy(data,MOD,"CPF"))

  MOD$tau <- CTMM$tau[2]
  MODELS <- c(MODELS,ctmm.select.isotropy(data,MOD,"OU"))

  MOD <- CTMM$tau <- NULL
  MODELS <- c(MODELS,ctmm.select.isotropy(data,MOD,"IID"))

  return(MODELS)  
}


ctmm.select.ouf <- function(data,CTMM)
{
  MODELS <- list()
  NAMES <- c("IID","OU","OUF")

  # step down in autocorrelation degree
  for(K in length(CTMM$tau):0) 
  {
    MOD <- CTMM

    # fix K
    if(K == 0) { MOD$tau <- NULL } else { MOD$tau <- MOD$tau[1:K] }

    MODELS <- c(MODELS,ctmm.select.isotropy(data,MOD,NAMES[K+1]))
  }

  return(MODELS)
}


ctmm.select.isotropy <- function(data,CTMM,NAME=NULL)
{
  MODELS <- list()
  MOD <- CTMM
  
  # name model
  if(MOD$isotropic)
  { name <- paste(NAME,"isotropic") }
  else
  { name <- paste(NAME,"anisotropic") } 
  
  # given model
  MODELS <- c(MODELS,ctmm.select.circle(data,MOD,name))
  
  # step down in isotropy
  if(CTMM$isotropic==FALSE)
  {
    MOD$isotropic <- TRUE
    MOD$sigma <- covm(MOD$sigma,isotropic=TRUE)
    name <- paste(NAME,"isotropic")
    MODELS <- c(MODELS,ctmm.select.circle(data,MOD,name))
  }
  
  return(MODELS)
}


ctmm.select.circle <- function(data,CTMM,NAME=NULL)
{
  MODELS <- list()
  MOD <- CTMM
  
  if(length(CTMM$tau)==0 || CTMM$tau==0)
  { CTMM$circle <- FALSE }
  
  # given model
  if(MOD$circle)
  { name <- paste(NAME,"circle") }
  else
  { name <- NAME }
  MODELS[[name]] <- MOD
  
  #simpler alternative
  if(MOD$circle)
  {
    MOD$circle <- FALSE
    MODELS[[NAME]] <- MOD
  }
  
  return(MODELS)
}


#################################################
# SLOW LIKELIHOOD FUNCTION
# not to be exposed to end users
# THIS NEEDS TO BE FIXED WITH ERRORS
ctmm.loglike.slow <- function(data,CTMM)
{  
  t <- data$t
  x <- data$x
  y <- data$y

  tau <- CTMM$tau
  sigma <- CTMM$sigma
  mu <- CTMM$mu
  isotropic <- CTMM$isotropic
  
  K <- length(tau)
  n <- length(t)
  
  # lag matrix
  C <- outer(t,t,"-")
  C <- abs(C)
  
  # now the correlation matrix
  if(K==1)
  { C <- exp(-C/tau) }
  else if(K==2)
  { C <- (tau[1]*exp(-C/tau[1])-tau[2]*exp(-C/tau[2]))/(tau[1]-tau[2]) }
  
  logdetC <- determinant(C,logarithm=TRUE)$modulus
  
  u <- rep(1,n) # stationary mean function
  
  # now the inverse correlation matrix times vectors
  C <- solve(C,cbind(x,y,u))
  
  Cx <- C[,1]
  Cy <- C[,2]
  w <- C[,3]

  W <- sum(w)
  
  if(is.null(mu)) { mu <- c(sum(Cx), sum(Cy))/W }
  
  # detrend mean
  x <- x - mu[1]
  y <- y - mu[2]

  Cx <- Cx - mu[1]*w
  Cy <- Cy - mu[2]*w
  
  if(is.null(sigma))
  {
    sigma <- x %*% Cy
    sigma <- rbind( c(x %*% Cx, sigma) , c( sigma, y %*% Cy) )/n
    sigma <- covm(sigma,isotropic)
  }

  if(length(sigma) == 1) { sigma <- sigma * diag(2)  }
  
  COV.mu <- methods::getDataPart(sigma)/W
  
  loglike <- -logdetC - (n/2)*log(det(2*pi*sigma)) - n

  CTMM <- ctmm(loglike=loglike,tau=tau,sigma=sigma,mu=mu,COV.mu=COV.mu,DOF.mu=W,info=attr(data,"info"))
    
  return(CTMM)
}


###################
# general parameter guessing function
###################
ctmm.guess <- function(data,variogram=NULL,CTMM=ctmm(),name="GUESS",interactive=TRUE)
{
  if(is.null(variogram)) { variogram = variogram(data) }
  
  mu <- c(mean(data$x),mean(data$y))
  
  z <- cbind(data$x - mu[1], data$y - mu[2])
  sigma <- t(z) %*% z
  
  # remove error from variability
  if(CTMM$error & !is.null(data$error)) { sigma <- sigma - sum(data$error) * diag(2) }
  
  n <- length(data$t)
  sigma <- sigma / (n-1)  

  CTMM$mu <- mu
  CTMM$sigma <- covm(sigma)
  
  # estimate circulation period
  if(CTMM$circle==1)
  {
    # velocities
    v <- cbind(diff(z[,1]),diff(z[,2])) / diff(data$t)
    # midpoint locations during velocity v
    z <- cbind(z[-1,1]+z[-n,1],z[-1,2]+z[-n,2])/2

    # average angular momentum
    L <- (z[,1]%*%v[,2] - z[,2]%*%v[,1]) / (n-1)
    
    circle <- L / mean(diag(sigma))
    circle <- 2*pi/circle
    
    CTMM$circle <- circle
  }
  
  variogram.fit(variogram,CTMM=CTMM,name=name,interactive=interactive)
}
