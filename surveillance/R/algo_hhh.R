###################################################
### chunk number 1: 
###################################################
# lag - which lag for observation-driven part?
#  y_i,t = lambda*y_i,t-lag  NOTE: lag=-1 means y_i,t+1
# lag.range =c(lag.neg, lag.pos) 
#  i.e. (1,0) for t-1,t (DEFAULT)
#       (2,2) for t-2,t-1,t,t+1,t+2
algo.hhh <- function(disProgObj, 
                     control=list(lambda=TRUE, 
                                  neighbours=FALSE, 
                                  linear=FALSE, 
                                  nseason=0,
                                  negbin=c("none", "single", "multiple"), 
                                  proportion=c("none", "single", "multiple"),
                                  lag.range=NULL),
                     thetastart=NULL, 
                     verbose=TRUE){
        
  #Convert sts objects
  if (class(disProgObj) == "sts") disProgObj <- sts2disProg(disProgObj)

  #set default values (if not provided in control)
  if(is.null(control[["linear",exact=TRUE]]))
    control$linear <- FALSE
    
  if(is.null(control[["nseason",exact=TRUE]]))
    control$nseason <- 0
    
  if(is.null(control[["neighbours",exact=TRUE]]))
    control$neighbours <- NA
    
  if(is.null(control[["negbin",exact=TRUE]]))
    control$negbin <- "none"
    
  if(is.null(control[["lambda",exact=TRUE]]))
    control$lambda <- 1

  if(is.null(control[["proportion",exact=TRUE]]))
    control$proportion <- "none"
    
  control$negbin <- match.arg(control$negbin, c("single","multiple","none"))
  control$proportion <- match.arg(control$proportion, c("single","multiple","none"))

  # convert logical values to numerical values, FALSE corresponds to NA
  # to allow for lag == 0
  if(is.logical(control[["lambda", exact=TRUE]])){
    control$lambda <- as.numeric(control$lambda)
    control$lambda[control$lambda==0] <- NA
  }
   if(is.logical(control[["neighbours", exact=TRUE]])){
    control$neighbours <- as.numeric(control$neighbours)
    control$neighbours[control$neighbours==0] <- NA
  }
  
  # determine range of observations y_i,t
  if(is.null(control[["lag.range",exact=TRUE]])){
    lags <- c(control$lambda, control$neighbours)
    control$lag.range <- c(max(c(lags,1),na.rm=TRUE), 
                           max(c(-lags,0), na.rm=TRUE))
  }

  # check if observed is a vector and convert to matrix if necessary
  if(is.vector(disProgObj$observed)) disProgObj$observed <- as.matrix(disProgObj$observed)
  n <- nrow(disProgObj$observed)
  nareas <- ncol(disProgObj$observed)

  #univariate
  if(nareas ==1){
    control$neighbours <- NA
    control$proportion <- "none"
    
    control$nseason <- control$nseason[1]
  }

  # model with (lambda, pi) ?
  if(control$proportion != "none"){
    control$neighbours <- NA
    # no lambda specified or lambda not specified for each area
    if(sum(!is.na(control$lambda)) == 0 | sum(!is.na(control$lambda))!= nareas)
      control$lambda <- 1
  }
  
  # check neighbourhood matrix if neighbours=TRUE or proportion!="none"
  if(sum(!is.na(control$neighbours))>0 | control$proportion != "none"){
    # is there a neighbourhood matrix?
    if(is.null(disProgObj$neighbourhood)) stop("No neighbourhood matrix is provided\n")
    if(any(is.na(disProgObj$neighbourhood))) stop("No correct neighbourhood matrix given\n")
  }

  #make "design" matrices
  designRes<- make.design(disProgObj=disProgObj, control=control)

  # check if there are neighbours
  if(designRes$dim$phi > 0){
    nOfNeighbours <- designRes$nOfNeighbours
    
    if((designRes$dim$phi==1) & (sum(nOfNeighbours)==0))
      stop("Specified model is not in line with neighbourhood matrix\n")
      
#    if((designRes$dim$phi==nareas) & (any(nOfNeighbours[!is.na(control$neighbours)]==0)))
    if((length(control$neighbours) == nareas) & (any(nOfNeighbours[!is.na(control$neighbours)]==0)))
      stop("Specified model is not in line with neighbourhood matrix\n")
      
  } else if(designRes$dim$proportion > 0){
    nOfNeighbours <- designRes$nOfNeighbours
    
    if((designRes$dim$proportion==1) & (sum(nOfNeighbours)==0))
      stop("Specified model is not in line with neighbourhood matrix\n")
      
    if((designRes$dim$proportion==nareas) & (any(nOfNeighbours==0)))
      stop("Specified model is not in line with neighbourhood matrix\n")
      
  }


  dimtheta <- designRes$dimTheta$dim
  dimLambda <- designRes$dimTheta$lambda
  dimPhi <- designRes$dimTheta$phi

  #starting values for optim 
  areastart <- log(colMeans(designRes$Y)/designRes$populationFrac[1,])
  
  if(!is.null(thetastart)){
    #check dimension of thetastart
    # must be either of length dimtheta or of length dimtheta-nareas
    if(all(length(thetastart) != c(dimtheta, dimtheta-nareas)) ){
      cat('thetastart must be of length', dimtheta, 'or ',dimtheta-nareas,'\n')
     return(NULL)
    }
    theta  <- thetastart
    if(length(theta) == dimtheta) areastart <- NULL
  } else {
    #set starting values for theta
    #lambda = log(0.5), phi = log(0.1), beta = gamma = delta = 0, psi = 1
    theta <- c(rep(log(0.5),designRes$dimTheta$lambda), rep(log(0.1),designRes$dimTheta$phi), 
               rep(0.5,designRes$dimTheta$proportion),
               rep(0, designRes$dimTheta$trend + designRes$dimTheta$season), rep(2,designRes$dimTheta$negbin) )
  }
  
  #starting values for intercepts
  if(!is.null(areastart)){
    if(dimLambda + dimPhi >0){
      #cat("theta",theta[1:(dimLambda + dimPhi)],"\n")
      Lambda <- getLambda(theta[1:(dimLambda + dimPhi)], designRes)
      expAlpha <- expAlpha.mm(Lambda,designRes$Y)
      expAlpha[expAlpha <=0] <- (colMeans(designRes$Y)/designRes$populationFrac[1,])[expAlpha <=0]
      areastart <- log(expAlpha)
      #areastart <- log(expAlpha.mm(Lambda,designRes$Y))
    }
  
    theta <- c(areastart,theta)
    #cat("initial values",theta,"\n")
  }

  #check if initial values are valid
  mu<-meanResponse(theta,designRes)$mean
  if(any(mu==0) | any(!is.finite(mu)))
    stop("invalid initial values\n")

  # maximize loglikelihood
  mycontrol <- list(fnscale=-1, type=3, maxit=1000)
  suppressWarnings(myoptim <- optim(theta, fn=loglikelihood, gr=gradient,
                                  control=mycontrol, method="BFGS", hessian=TRUE, designRes=designRes))

  if(myoptim$convergence==0){
    convergence <- TRUE
  } else {
      if(verbose)
        cat("Algorithm has NOT converged. \n")
      res <- list(convergence=FALSE)
      class(res) <- "ah"
      return(res)
  }
  
  loglik <- myoptim$value
  
  if(loglik==0){
    if(verbose){
      cat('loglikelihood = 0\n')
      cat('Results are not reliable! Try different starting values. \n')
    }
    res <- list(convergence=FALSE)
    class(res) <- "ah"
    return(res)
  }

  thetahat <- myoptim$par
  fisher <- -myoptim$hessian

  # fitted values
  fitted <- meanResponse(thetahat,designRes)$mean     

  #psi, lambda and phi are on log-scale
  #-> transformation of estimates, standard errors and fisher (using delta rule)
  #labels for results
  D <- jacobian(thetahat, designRes)$D
  thetahat <- jacobian(thetahat, designRes)$theta

  #Approximation to the inverted fisher info matrix
#  cov <- try(D %*% solve(fisher) %*% D, silent=TRUE)
  cov <- try(D %*% solve(fisher) %*% t(D), silent=TRUE)

  #fisher info is singular
  if(class(cov) == "try-error"){
    if(verbose){
      cat("Fisher info singular \t loglik=",loglik," \n")
      cat("theta",round(thetahat,2),"\n")
      cat('Results are not reliable! Try different starting values. \n')
    }
    res <- list(convergence=FALSE)
    class(res) <- "ah"
    return(res)
  }

  if(any(!is.finite(diag(cov))) | any(diag(cov)<0)){
    if(verbose){
      cat("infinite or negative cov\t loglik=",loglik,"\n")
      cat("theta",round(thetahat,2),"\n")
      cat('Results are not reliable! Try different starting values. \n')
    }
    res <- list(convergence=FALSE)
    class(res) <- "ah"
    return(res)
  }

  se <- sqrt(diag(cov))

  if(convergence & verbose)
    cat("Algorithm claims to have converged \n")
  
  result <- list(coefficients=thetahat, se=se, cov=cov, call=match.call(),
                 loglikelihood=loglik, convergence=convergence,
                 fitted.values=fitted, control=control,disProgObj=disProgObj, 
                 lag=designRes$lag, nObs=designRes$nObs)
  
  class(result) <- "ah"
  return(result)
}



###################################################
### chunk number 2: 
###################################################
algo.hhh.grid <- function(disProgObj, control=list(lambda=TRUE,neighbours=FALSE, 
               linear=FALSE, nseason=0, 
               negbin=c("none", "single", "multiple"), 
               proportion=c("none", "single", "multiple"),lag.range=NULL), 
               thetastartMatrix, maxTime=1800, verbose=FALSE){

  #convert disProgObj if necessary
  if (class(disProgObj) == "sts") disProgObj <- sts2disProg(disProgObj)

  #set default values (if not provided in control)
  if(is.null(control[["linear",exact=TRUE]]))
    control$linear <- FALSE
    
  if(is.null(control[["nseason",exact=TRUE]]))
    control$nseason <- 0
    
  if(is.null(control[["neighbours",exact=TRUE]]))
    control$neighbours <- NA
    
  if(is.null(control[["negbin",exact=TRUE]]))
    control$negbin <- "none"
    
  if(is.null(control[["lambda",exact=TRUE]]))
    control$lambda <- 1

  if(is.null(control[["proportion",exact=TRUE]]))
    control$proportion <- "none"
    
  control$negbin <- match.arg(control$negbin, c("single","multiple","none"))
  control$proportion <- match.arg(control$proportion, c("single","multiple","none"))
  
  # convert logical values to numerical values, FALSE corresponds to NA
  # to allow for lag == 0
  if(is.logical(control[["lambda", exact=TRUE]])){
    control$lambda <- as.numeric(control$lambda)
    control$lambda[control$lambda==0] <- NA
  }
   if(is.logical(control[["neighbours", exact=TRUE]])){
    control$neighbours <- as.numeric(control$neighbours)
    control$neighbours[control$neighbours==0] <- NA
  }
  
  # determine range of observations y_i,t
  if(is.null(control[["lag.range",exact=TRUE]])){
    lags <- c(control$lambda, control$neighbours)
    control$lag.range <- c(max(c(lags,1),na.rm=TRUE), 
                           max(c(-lags,0), na.rm=TRUE))
  } 
  n <- nrow(disProgObj$observed)
  nareas <- ncol(disProgObj$observed)

  # check parameter specification for season
  #univariate
  if(nareas ==1){
    control$neighbours <- NA
    control$proportion <- "none"
    
    control$nseason <- control$nseason[1]
  }

  # model with (lambda, pi) ?
  if(control$proportion != "none"){
    control$neighbours <- NA
    # no lambda specified or lambda not specified for each area
    if(sum(!is.na(control$lambda)) == 0 | sum(!is.na(control$lambda))!= nareas)
      control$lambda <- 1
  }
  
  # check neighbourhood matrix if neighbours=TRUE or proportion!="none"
  if(sum(!is.na(control$neighbours))>0 | control$proportion != "none"){
    if(any(is.na(disProgObj$neighbourhood)))
    stop("No correct neighbourhood matrix given\n")
  }

  designRes<- make.design(disProgObj=disProgObj, control=control)

  # check if there are neighbours
  if(designRes$dim$phi > 0){
    nOfNeighbours <- designRes$nOfNeighbours
    
    if((designRes$dim$phi==1) & (sum(nOfNeighbours)==0))
      stop("Specified model is not in line with neighbourhood matrix\n")
      
#    if((designRes$dim$phi==nareas) & (any(nOfNeighbours[!is.na(control$neighbours)]==0)))
    if((length(control$neighbours) == nareas) & (any(nOfNeighbours[!is.na(control$neighbours)]==0)))
      stop("Specified model is not in line with neighbourhood matrix\n")
      
  } else if(designRes$dim$proportion > 0){
    nOfNeighbours <- designRes$nOfNeighbours
    
    if((designRes$dim$proportion==1) & (sum(nOfNeighbours)==0))
      stop("Specified model is not in line with neighbourhood matrix\n")
      
    if((designRes$dim$proportion==nareas) & (any(nOfNeighbours==0)))
      stop("Specified model is not in line with neighbourhood matrix\n")
      
  }
  
  dimthetaStart <- designRes$dimTheta$dim -nareas
  
  if(dimthetaStart == 0){ #only intercepts, grid search not necessary                    
    return(algo.hhh(disProgObj=disProgObj,control=control))
  }
  
  #check dimension of thetastartMatrix
  if(!is.matrix(thetastartMatrix)){
    cat('thetastart must be a matrix with', designRes$dimTheta$dim, 'or ', dimthetaStart, 'columns\n')
    return(NULL)
  }
  if(all(ncol(thetastartMatrix) != c(designRes$dimTheta$dim, dimthetaStart))){
    cat('thetastart must be a matrix with', designRes$dimTheta$dim, 'or ', dimthetaStart,'columns\n')
    return(NULL)
  }
  
  #try multiple starting values and return the result with highest likelihood
  #stop search once time limit is exceeded
  i<-0
  nOfIter <- nrow(thetastartMatrix)
  gridUsed <- nOfIter
  
  if(verbose) cat('The size of grid is', nOfIter, '\n')
  
  bestLoglik <- list(loglikelihood = -1e99)
  allLoglik <- matrix(NA,nrow=nOfIter,ncol=designRes$dimTheta$dim+1)
  time <- maxTime
  
  while((time > 0) & (i < nOfIter)){
    i <- i+1
    #run algo.hhh with the i-th row of thetastartMatrix as initial values
    time.i <- system.time(res<-try(algo.hhh(disProgObj=disProgObj,control=control,thetastart=thetastartMatrix[i,],verbose=verbose),silent=!verbose))[3]
    #how much time is left now
    time <- time - time.i
    
    #print progress information
    if(verbose){
      if(class(res)== "try-error")
        print(c(niter=i,timeLeft=time,loglik=NULL))
      else print(c(niter=i,timeLeft=time,loglik=res$loglikelihood))
    }

    
    #don't consider "useless" results for the search of the best loglikelihood
    if(class(res)!= "try-error" && res$convergence){ 
      #save loglik
      allLoglik[i,] <- c(res$loglikelihood,coef(res))
      #keep it as bestLoglik if loglikelihood improved  
      if(res$loglikelihood > bestLoglik$loglikelihood){
        bestLoglik <- res
      }
    }
  }
  
  if(time < 0){
    if(verbose) cat('Time limit exceeded, grid search stopped after', i, 'iterations. \n')
    allLoglik <- as.matrix(allLoglik[1:i,])
    gridUsed <- i
  }
  
  timeUsed <- ifelse(time>0, maxTime-time,maxTime+abs(time))
  
  #algo.hhh did not converge or produced useless results for all starting values, 
  #i.e. there is no result
  if(is.null(coef(bestLoglik))) {
    #convergence <- FALSE
    #cat('Algorithms did not converge, please try different starting values! \n')
    bestLoglik <- list(loglikelihood=NULL,convergence=FALSE)
  } else{
  #give names to all Loglik-matrix
  colnames(allLoglik) <- c("loglik",names(coef(bestLoglik)))
  }
  
  result <- list(best = bestLoglik, allLoglik = allLoglik,gridSize=nOfIter,gridUsed=gridUsed, time=timeUsed,convergence=bestLoglik$convergence)
  class(result) <- "ahg"
  return(result) 

}



###################################################
### chunk number 3: 
###################################################
create.grid <- function(disProgObj, control, params = list(epidemic = c(0.1, 0.9, 5),
                 endemic=c(-0.5,0.5,3), negbin = c(0.3, 12, 10))) {
  
  #convert S4 sts to S3 if necessary
  if (class(disProgObj) == "sts") disProgObj <- sts2disProg(disProgObj)

  designRes <- make.design(disProgObj, control)
  control <- designRes$control
  dimParams <- designRes$dimTheta
  
  dimLambda <- dimParams$lambda
  dimPhi <- dimParams$phi
  dimProp <- dimParams$proportion
  
  dimEndemic <- dimParams$trend + dimParams$season
  
  dimNegbin <- dimParams$negbin
  
  # check if initial values are provided
  if((dimLambda +dimPhi > 0) & is.null(params$epidemic))
    stop("Please provide initial values for the epidemic component \n")
  if((dimEndemic > 0) & is.null(params$endemic))
    stop("Please provide initial values for the endemic component \n")
  if((dimNegbin > 0) & is.null(params$negbin))
    stop("Please provide initial values for the dispersion parameter psi \n")
  
  # check if initial values are specified correctly
  if(!is.null(params$epidemic)){
   if( params$epidemic[3]%%1 !=0 | params$epidemic[3]<1 | sign(params$epidemic[3])== -1)
    stop("Last component of params$epidemic must be a positive integer\n")
  }
  if(!is.null(params$endemic)){ 
   if( params$endemic[3]%%1 !=0 | params$endemic[3]<1 | sign(params$endemic[3])== -1)
    stop("Last component of params$endemic must be a positive integer\n")
  }
  if(!is.null(params$negbin)){ 
   if( params$negbin[3]%%1 !=0 | params$negbin[3]<1 | sign(params$negbin[3])== -1)
    stop("Last component of params$negbin must be a positive integer\n")
  }
    
  grid <- list()

  if(dimNegbin >0){
    psi <- seq(params$negbin[1], params$negbin[2], length = params$negbin[3])
    if(any(psi<=0))
      stop("Initial values for psi must be positive\n")
    log.psi <- log(psi[psi >0])
    grid$psi <- log.psi
  }
  
  if(dimLambda >0){
    epidemic <- seq(params$epidemic[1], params$epidemic[2], length = params$epidemic[3])
    if(any(epidemic<=0))
      stop("Iinitial values for lambda must be positive\n")
    log.lambda <-  log(epidemic[epidemic >0])
    grid$lambda <- log.lambda
  }
  
  if(dimPhi >0){
    epidemic <- seq(params$epidemic[1], params$epidemic[2], length = params$epidemic[3])
    if(any(epidemic<=0))
      stop("Initial values for phi must be positive\n")
    log.lambda <-  log(epidemic[epidemic >0])
    grid$phi <- log.lambda
  }
  
  if(dimProp >0){
    if(any(epidemic<=0 | epidemic >=1))
      stop("initial values for pi must be in (0,1)\n")
    logit.prop <- log(epidemic[epidemic > 0 & epidemic < 1]) - log(1-epidemic[epidemic > 0 & epidemic < 1])
    grid$prop <- logit.prop
  }
  
  if(dimEndemic >0){
    endemic <- seq(params$endemic[1], params$endemic[2], length = params$endemic[3])
    grid$endemic <- endemic
  }
    
  grid <- expand.grid(grid)
  
  grid <- as.matrix( grid[c(rep("lambda",dimLambda), rep("phi",dimPhi), rep("prop",dimProp),
                          rep("endemic",dimEndemic), rep("psi",dimNegbin))] )
                          
  gridSize <- nrow(grid)
  cat("Matrix with ",gridSize, " initial values\n")
  
  return(grid)
}




###################################################
### chunk number 4: 
###################################################
loglikelihood <- function(theta, designRes){
  
  control <- designRes$control
  Y <- designRes$Y

  mean <- meanResponse(theta=theta, designRes=designRes)$mean
  
  dimNegbin <- designRes$dimTheta$negbin
  dimTheta <- designRes$dimTheta$dim

  #loglikelihood poisson
  if(dimNegbin==0){
    result <- colSums(dpois(Y, lambda=mean, log=TRUE))
  
  } else if(dimNegbin==1){
    #loglikelihood negbin
    
    #ensure psi (on last position in vector theta) ist positive
    psi <- exp(theta[dimTheta])

    result <- colSums(dnbinom(Y, size=psi, mu=mean, log=TRUE))
    
  } else if(dimNegbin>1){
    #loglikelihood negbin, multiple dispersion params
    
    #ensure psi (on last positions) is positive
    psi <- exp(theta[(dimTheta-dimNegbin+1):dimTheta])
    psi <- matrix(psi,ncol=dimNegbin, nrow=nrow(Y), byrow=TRUE)

    result <- colSums(dnbinom(Y, size=psi, mu=mean, log=TRUE))
  }
  res <- sum(result)
  attr(res, "colsums") <- result
  return(res)
}



###################################################
### chunk number 5: 
###################################################
meanResponse <- function(theta, designRes){

  # unpack design matrices
  Y <- designRes$Y  
  nareas <- ncol(Y)
  n <- nrow(Y)
  
  X.trendSeason <- designRes$X.trendSeason
  Ym1 <- designRes$Ym1
  
  Ym1.neighbours <- designRes$Ym1.neighbours
  nhood <- designRes$disProgObj$neighbourhood
  nOfNeighbours <- designRes$nOfNeighbours
  
  pop <- designRes$populationFrac

  #check dimension of theta
  if(designRes$dimTheta$dim != length(theta)){
    cat('theta must be of length',designRes$dimTheta$dim,'\n')
    return(NULL)
  }

  #unpack parameters and ensure lambda and phi are positive
  params <- unpackParams(theta,designRes)


  ###################################################################
  ## calculation of mean

  #autoregressive part
  
  # model with lambda and phi ?
  if(designRes$control$proportion == "none"){

    #auto=0 if lambda and phi are not used in model
    lambda <- params$lambda
    phi <- params$phi
    # no autoregression
    if(is.null(lambda))
      auto.lambda<- 0
    else {
      # multiple lambda
      if(length(designRes$control$lambda)==nareas){
        # create vector lambda with elements 0 if control$lambda=FALSE
        lambda <- rep(0,nareas)
        lambda[!is.na(designRes$control$lambda)] <- params$lambda
      }
      auto.lambda <- Ym1*matrix(lambda,ncol=nareas,nrow=nrow(Y), byrow=TRUE)
    }

    if(is.null(phi))
       auto.phi <- 0
    else {
      # multiple phi
      if(length(designRes$control$neighbours)==nareas){
        # create vector phi with elements 0 if control$neighbours=FALSE
        phi <- rep(0,nareas)
        phi[!is.na(designRes$control$neighbours)] <- params$phi
      }
      
      auto.phi <- Ym1.neighbours*matrix(phi,ncol=nareas,nrow=nrow(Y), byrow=TRUE)
    }

    auto<- auto.lambda + auto.phi
  } else {
    #################################################
    ## model with lambda and proportion pi
    #################################################
    
    # helper function
    weightedSumEpidemic <- function(prop,lambda){
    
      # ensure region id is not included
      diag(nhood) <- 0
      # compute sum_j~i {pi_ji * Y_j,t-1} for unit id
      # where pi_ji = pi_i          for j=i
      #       pi_ji =(1-pi_j)/|k~j| for j~i
      
      one <- function(id){
        # nOfNeighbours is number of Neigbours for each unit id=1,..,m  i.e. |k~id| 
        nOfNeighbours[id]<-1
        pi.ij <- matrix(lambda*(1-prop)/nOfNeighbours,ncol=length(prop),nrow=nrow(Ym1),byrow=TRUE)
        # select pi_ij with j~i
        piYm1 <- as.matrix((Ym1*pi.ij)[,nhood[,id]>0])
        rowSums(piYm1)
      }
      sapply(1:length(prop),one)
    }

    lambda <- matrix(params$lambda,ncol=nareas,nrow=n,byrow=TRUE)
    
    if(designRes$control$proportion == "single")
      prop <- rep(params$pi,nareas)
    else prop <- params$pi
    
    # lambda*pi_ji*Y_j,t-1
    auto.phi <- weightedSumEpidemic(prop=prop,lambda=lambda[1,])

    auto.lambda <-  Ym1*lambda*matrix(prop,ncol=nareas,nrow=n,byrow=TRUE)
    
    auto <- auto.lambda+auto.phi
  }
  ################
  
  #trend and seasonal components
  nSeason <- designRes$control$nseason
  dimSeason <- designRes$dimTheta$season
  dimTrend <- designRes$dimTheta$trend
  
  # trend
  if(dimTrend >0){
    if(length(designRes$control$linear) == 1)
      beta <- rep(params$beta,nareas)
    else {
      beta <- rep(0,nareas)
      beta[designRes$control$linear] <-params$beta
    }
    
    predTime <- as.matrix(X.trendSeason[,1])%*%beta
  } else predTime <- 0
  
  # season
  if( dimSeason >0){
    # discard design matrix for trend
    X.Season <- X.trendSeason[,(1+ (dimTrend>0) ):ncol(X.trendSeason)]
    maxSeason <- max(nSeason)

    #construct a suitable matrix for seasonal parameters gamma_i
    # same number of Fourier frequencies S for all areas i:
    if(length(nSeason)==1){
      gammaMat <- matrix(params$gamma,ncol=nareas,nrow=2*maxSeason,byrow=FALSE)
      
    } else if(length(nSeason)==nareas){
    # different number of frequencies S_i for each area
      gammaMat <- matrix(0,ncol=nareas,nrow=2*maxSeason)
      index <- rep(1:nareas,2*nSeason)
      
      for(i in 1:nareas)
        gammaMat[0:(2*nSeason[i]),i] <- params$gamma[index==i]

    } else stop("nseason must be a vector of length 1 or",nareas,"\n")
    
    predSeason <- X.Season%*%gammaMat

      
  } else  predSeason <- 0 

  #intercepts for areas
  #matrix with columns (alpha_1,...,alpha_nareas)
  predarea <- matrix(params$alpha, byrow=TRUE, ncol=nareas, nrow=nrow(Y))
  
  #endemic part 
  endemic <- pop*exp(predarea+predTime+predSeason)
  #results
  mean <- auto + endemic

  #Done
  return(list(mean=mean,epidemic=auto,endemic=endemic,epi.own=auto.lambda,epi.neighbours=auto.phi))
}



###################################################
### chunk number 6: 
###################################################
make.design <- function(disProgObj, control=list(lambda=TRUE, neighbours=FALSE, 
        linear=FALSE, nseason=0,
         negbin=c("none", "single", "multiple"), 
         proportion=c("none", "single", "multiple"),
         lag.range=NULL) ){

  #Convert sts objects
  if (class(disProgObj) == "sts") disProgObj <- sts2disProg(disProgObj)

  #set default values (if not provided in control)
  if(is.null(control[["linear",exact=TRUE]]))
    control$linear <- FALSE
    
  if(is.null(control[["nseason",exact=TRUE]]))
    control$nseason <- 0
  
  if(is.null(control[["neighbours",exact=TRUE]]))
    control$neighbours <- NA
    
  if(is.null(control[["negbin",exact=TRUE]]))
    control$negbin <- "none"
    
  if(is.null(control[["lambda",exact=TRUE]]))
    control$lambda <- 1
    
  if(is.null(control[["proportion",exact=TRUE]]))
    control$proportion <- "none"
    
  control$proportion <- match.arg(control$proportion, c("single","multiple","none"))
  control$negbin <- match.arg(control$negbin, c("single","multiple","none"))
  
  # convert logical values to numerical values, FALSE corresponds to NA
  # to allow for lag == 0
  if(is.logical(control[["lambda", exact=TRUE]])){
    control$lambda <- as.numeric(control$lambda)
    control$lambda[control$lambda==0] <- NA
  }
   if(is.logical(control[["neighbours", exact=TRUE]])){
    control$neighbours <- as.numeric(control$neighbours)
    control$neighbours[control$neighbours==0] <- NA
  }
  
  # determine range of observations y_i,t
  if(is.null(control[["lag.range",exact=TRUE]])){
    lags <- c(control$lambda, control$neighbours)
    control$lag.range <- c(max(c(lags,1),na.rm=TRUE), 
                           max(c(-lags,0), na.rm=TRUE))
  }
    
  data <- disProgObj$observed
  n <- nrow(data)
  nareas <- ncol(data)
  
  # check parameters
  if(length(control$lambda)>1 & length(control$lambda)!=nareas)
    stop("parameter lambda is not specified correctly\n")
    
  if(length(control$neighbours)>1 & length(control$neighbours)!=nareas)
    stop("parameter phi is not specified correctly\n")
    
  if(length(control$linear)>1 & length(control$linear)!=nareas)
    stop("parameter beta is not specified correctly\n")
  
  #univariate
  if(nareas ==1){
    control$neighbours <- NA
    control$proportion <- "none"
    
    control$nseason <- control$nseason[1]
  }

  # maximum number of seasonal Fourier frequencies
  maxSeason <- max(control$nseason)

  # model with (lambda, pi) ?
  if(control$proportion != "none"){
    control$neighbours <- NA
    # no lambda specified or lambda is not specified for each area
    if(sum(!is.na(control$lambda)) == 0 |sum(!is.na(control$lambda)) !=nareas)
      control$lambda <- 1
  }

  dimLambda <- sum(!is.na(control$lambda))

  dimPhi <- sum(!is.na(control$neighbours))

  dimProportion <- switch(control$proportion ,
                       "single" = 1,
                       "multiple"= nareas,
                       "none" = 0)

  dimTrend <-  sum(control$linear)
  
  dimSeason <- sum(2*control$nseason)
  
  dimIntercept <- nareas

  dimNegbin <- switch(control$negbin, 
                       "single" = 1,
                       "multiple"= nareas,
                       "none" = 0)

  #theta = (alpha_i, lambda, phi (or pi), beta,gamma_i,delta_i,..., psi)
  dim <- dimLambda+dimPhi+dimTrend+dimSeason+dimIntercept+dimProportion+dimNegbin
  dimTheta <- list(lambda=dimLambda, phi=dimPhi, trend=dimTrend, season=dimSeason,
                   intercept=dimIntercept, proportion=dimProportion, negbin=dimNegbin ,dim=dim)

  ####################################################################
  # arrange response as matrix
  #Y, Ym1, Ym1.neighbours and population are (nOfobs)x(nOfareas) matrices
  #where nOfareas is the number of areas/units and 
  # nOfobs is determined by control$lag.range with default nOfObs=n-1 
  
  # Thus, lag.range can be used to ensure that models with different lags
  # are based on the same observations.
  t.min <- 1+control$lag.range[1]
  t.max <- n-control$lag.range[2]

  Y <- matrix(data[t.min:t.max,],nrow=length(t.min:t.max),ncol=nareas)
  
  # population sizes n_{i,t}
  if(is.null(disProgObj$populationFrac)){
	population <- matrix(1, nrow=length(t.min:t.max),ncol=nareas)
  } else {
	population <- matrix(disProgObj$populationFrac[t.min:t.max,],nrow=length(t.min:t.max),ncol=nareas)
  }

  # observed counts at time point t-lag
  # NOTE: the same lag (the maximum lag) is used for all areas
  if(dimLambda >0){
    lag.lambda <- control$lambda[which.max(abs(control$lambda))]
    Ym1 <- matrix(data[(t.min:t.max)-lag.lambda,],nrow=length(t.min:t.max),ncol=nareas)
  } else {
    lag.lambda<- NA
    Ym1 <- matrix(0,nrow=length(t.min:t.max),ncol=nareas)
  }

  Ym1.neighbours <- matrix(0,nrow=length(t.min:t.max),ncol=nareas)
  nOfNeighbours <- 0

  # now matrix for neighbours
  if(dimPhi>0){
    lag.phi <- control$neighbours[which.max(abs(control$neighbours))]
    Ym1.neighbours <- weightedSumNeighbours(disProgObj)$neighbours[(t.min:t.max)-lag.phi,]
    nOfNeighbours <- weightedSumNeighbours(disProgObj)$nOfNeighbours
#     Ym1.neighbours <- sumNeighbours(disProgObj)[-n,]
 } else lag.phi <- NA
 
  if(dimProportion >0){
    Ym1.neighbours <- weightedSumNeighbours(disProgObj)$neighbours[(t.min:t.max)-lag.lambda,] #not really needed
    nOfNeighbours <- weightedSumNeighbours(disProgObj)$nOfNeighbours
  }
  
  ####################################################################
  # now define design matrix (for trend and seasonality) for each time point

  #t<- disProgObj$week[t.min:t.max]
  # if no $week is given 
  if(is.null(disProgObj$week)){
	t <- (t.min:t.max)-1
  } else {
	t<- disProgObj$week[(t.min:t.max)-1]
  }
  #t <- t - mean(t)

  form<-function(mod=ifelse(dimTrend == 0,"~-1","~-1+t"),
                 S=maxSeason, period=disProgObj$freq){
    if(S>0){
      for(i in 1:S){
        mod <- paste(mod,"+sin(",2*i,"*pi*t/",period,")+cos(",2*i,"*pi*t/",period,")",sep="")
      }
    }
    return(as.formula(mod))
  }

  if(dimTrend +dimSeason >0)
    X.trendSeason<-model.matrix(form(),data.frame(t=t))
  else X.trendSeason <-NULL


  result <- list("Y"=Y, "Ym1"=Ym1, "Ym1.neighbours"=Ym1.neighbours,"nOfNeighbours"=nOfNeighbours,
                 "X.trendSeason"=X.trendSeason,
                 "populationFrac"=population,  "dimTheta"=dimTheta,
                 "control"=control,"disProgObj"=disProgObj, "lag"=c(lag.lambda,lag.phi),"nObs"=prod(dim(Ym1)))

  return(result)
}



###################################################
### chunk number 7: 
###################################################
print.ah <- function(x,digits = max(3, getOption("digits") - 3), amplitudeShift=TRUE,reparamPsi=TRUE,...){
  if(!x$convergence)
    cat('Results are not reliable! Try different starting values. \n')
  else {
    if(!is.null(x$call)){
      cat("Call: \n")
      print(x$call)
    }
        
    cat('\nEstimated parameters and standard errors: \n\n')
    coefs <- coefficients(x, se=TRUE, amplitudeShift=amplitudeShift,reparamPsi=reparamPsi)
    
    print(round(cbind("Estimates"=coefs[,"Estimates"],
                 "Std.Error"=coefs[,"Std. Error"]),digits=digits),print.gap=2)

    cat('\nlog-likelihood:   ',round(x$loglik,digits=digits-2),'\n')  
    cat('AIC:              ',round(AIC(x),digits=digits-2),'\n')
    cat('BIC:              ',round(AIC(x,k=log(x$nObs)),digits=digits-2),'\n\n')
    
    if(!is.na(x$lag[1])) cat('lag used for lambda:      ',x$lag[1],'\n')
    if(!is.na(x$lag[2])) cat('lag used for phi:         ',x$lag[2] ,'\n')
    cat('number of observations:   ',x$nObs,'\n\n')
  }

}

print.ahg <- function (x, digits = max(3, getOption("digits") - 3), amplitudeShift=TRUE,reparamPsi=TRUE, ...){
    cat("\nsize of grid: ", x$gridSize, "\n")
    if (x$gridSize != x$gridUsed)
        cat("grid search stopped after", x$gridUsed, "iterations \n")
    cat("convergences: ",sum(!is.na(x$all[,1])),"\n")
    cat("time needed (in seconds)",x$time,"\n\n")
    if (!x$convergence)
        cat("\nAlgorithms did not converge, please try different starting values! \n")
    else {
      x$best$call <- NULL
      cat("values of log-likelihood:")
      print(table(round(x$all[,1],0)))
#      cat("\n")
      print.ah(x$best, digits = digits, amplitudeShift=amplitudeShift,reparamPsi=reparamPsi)
    }
}



###################################################
### chunk number 8: 
###################################################
#################################
# obtain predictions from the fitted algo.hhh model
#
# params:
#  object - a fitted object of class "ah" 
#  newdata - optionally, a disProgObject with which to predict; 
#            if omitted, the fitted mean is returned. 
#  type - the type of prediction required. The default is on the scale of the response 
#         variable (endemic and epidemic part) 
#         the alternative "endemic" returns only the endemic part (i.e. n_it * \nu_it)  
################################
predict.ah <- function(object,newdata=NULL,type=c("response","endemic","epi.own","epi.neighbours"),...){
  type <- match.arg(type,c("response","endemic","epi.own","epi.neighbours"))
  control <- object$control

  if(is.null(newdata))
    newdata <- object$disProgObj
  if(!inherits(newdata, "disProg"))
    stop("data must be an object of class disProg\n")

  coefs <- coefficients(object)

  design <- make.design(newdata,control=control)
  
  # in meanResponse the params lambda, phi are "exp()'ed"
  # log() them  to obtain the correct predictions
  if(sum(!is.na(control$lambda)) >0 | sum(!is.na(control$neighbours)) >0){
  
    indexL <- design$dimTheta$intercept+1
    indexU <- indexL +design$dimTheta$lambda +design$dimTheta$phi -1
    
    coefs[indexL:indexU] <- log(coefs[indexL:indexU])
    #cat("lambda,phi: indexL",indexL,"indexU",indexU,"\n")
    
    # pi is on logit-scale
    if(control$proportion != "none"){
      indexL <- design$dimTheta$intercept+design$dimTheta$lambda+1
      indexU <- indexL +design$dimTheta$proportion -1
      #cat("indexL",indexL,"indexU",indexU,"\n")
      coefs[indexL:indexU] <- log(coefs[indexL:indexU]/(1-coefs[indexL:indexU]))
    }

  }
  
  predicted <- meanResponse(coefs,design)

  if(type=="response")
    return(predicted$mean)
  else if(type=="endemic") return(predicted$endemic)
  else if(type=="epi.own") return(predicted$epi.own)
  else if(type=="epi.neighbours") return(predicted$epi.neighbours)
  
}

predict.ahg <- function(object, newdata=NULL, type=c("response","endemic","epi.own","epi.neighbours"),...){
  predict(object$best,newdata=newdata,type=type)
}




###################################################
### chunk number 9: 
###################################################
##########################
## residuals
##################
residuals.ah <- function (object, type = c("deviance", "pearson"), ...){
  type <- match.arg(type, c("deviance", "pearson"))
  
  # fitted values
  mean<- object$fitted.values
  #discard 1st observation (to obtain same dimension as mean)
  y <- as.matrix(object$disProgObj$observed[-1,])
    
  # poisson or negbin model
  if(object$control$negbin!="none"){
    coefs <- coefficients(object)
    psi <- matrix(coefs[grep("psi",names(coefs))],ncol=ncol(y),nrow=nrow(y),byrow=TRUE)
    
    distr <- function(mu){ dnbinom(y, mu=mu, size=psi, log=TRUE) }
    variance <- mean*(1+mean/psi)
  } else {
    distr <- function(mu){ dpois(y, lambda=mu,log=TRUE) }
    variance <- mean
  }

  res <- switch(type, 
                    deviance = sign(y-mean)*sqrt(2*(distr(y)-distr(mean))),
                    pearson = (y-mean)/sqrt(variance)
                )
  
  return(res)
}

residuals.ahg <- function(object, type = c("deviance", "pearson"), ...){
  residuals.ah(object$best,type=type)
}



###################################################
### chunk number 10: 
###################################################
############################################
# extract estimates and standard errors (se=TRUE)
# if amplitudeShift=TRUE, the seasonal params are transformed
# if reparamPsi=TRUE, the overdispersion param psi is transformed to 1/psi
#
############################################
coef.ah <- function(object,se=FALSE, amplitudeShift=FALSE, reparamPsi=FALSE,...){
  coefs <- object$coefficients
  stdErr <- object$se

  if(amplitudeShift & max(object$control$nseason)>0){
    #extract sin, cos coefficients
    index <- grep(" pi ",names(coefs))
    sinCos.names <- names(coefs)[index]
    # change labels
    names(coefs)[index] <- paste(c("A","s"),substr(sinCos.names,4,100),sep="")
    
    #transform sin, cos coefficients
    coefs[index] <- sinCos2amplitudeShift(coefs[index])
    # se's using Delta rule
    D <- diag(1,length(coefs))
    D[index,index]<- jacobianAmplitudeShift(coefs[index])
    cov <- D %*% object$cov %*% t(D)
    stdErr <- sqrt(diag(cov))
  }
  if(reparamPsi & object$control$negbin!="none"){
    #extract psi coefficients
    index <- grep("psi",names(coefs))
    psi.names <- names(coefs)[index]
    # change labels
    names(coefs)[index] <- paste("1/",psi.names,sep="")
    
    #transform psi coefficients
    coefs[index] <- 1/coefs[index]
    # se's using Delta rule: se[h(psi)] = se[psi] * |h'(psi)|
    # h = 1/psi, h' = -1/psi^2
    D <- diag(coefs[index]^2,length(index))
    stdErr[index] <- sqrt(diag(D %*% object$cov[index,index] %*% t(D)))
  }
  if(se)
    return(cbind("Estimates"=coefs,"Std. Error"=stdErr))
  else
    return(coefs)
}

coef.ahg <- function(object,se=FALSE, amplitudeShift=FALSE, reparamPsi=FALSE,...){
  return(coef(object$best,se=se, amplitudeShift=amplitudeShift,reparamPsi=reparamPsi))
}



###################################################
### chunk number 11: 
###################################################
## convert between sin/cos and amplitude/shift formulation
###################################################
# y = gamma*sin(omega*t)+delta*cos(omega*t)
#   =  A*sin(omega*t + phi)
# with Amplitude A= sqrt(gamma^2+delta^2)
# and shift phi= arctan(delta/gamma)
#################################################
sinCos2amplitudeShift <- function(params){
  # number of sin+cos terms
  lengthParams <- length(params)
  if(lengthParams %% 1 !=0)
    stop("wrong number of params")
  index.sin <- seq(1,lengthParams,by=2)
  
  one <- function(i=1,params){
    coef.sin <- params[i]
    coef.cos <- params[i+1]

    amplitude <- sqrt(coef.cos^2+coef.sin^2)
    shift <- atan2(coef.cos, coef.sin)
    return(c(amplitude,shift))
  }
  return(c(sapply(index.sin,one,params=params)))
}

amplitudeShift2sinCos <- function(params){
    lengthParams <- length(params)
    if (lengthParams%%1 != 0) 
        stop("wrong number of params")
    index.A <- seq(1, lengthParams, by = 2)
    one <- function(i = 1, params) {
        coef.A <- params[i]
        coef.shift <- params[i + 1]
        coef.cos <- -coef.A*tan(coef.shift)/sqrt(1+tan(coef.shift)^2)
        coef.sin <- -coef.A/sqrt(1+tan(coef.shift)^2)
        return(c(coef.sin,coef.cos))
    }
    return(c(sapply(index.A, one, params = params)))

}

##############################################
# y = gamma*sin(omega*t)+delta*cos(omega*t)
# g(gamma,delta) = [sqrt(gamma^2+delta^2), arctan(delta/gamma) ]'
# compute jacobian (dg_i(x)/dx_j)_ij
#############################################
jacobianAmplitudeShift <- function(params){
  # number of sin+cos terms
  lengthParams <- length(params)
  if(lengthParams %% 1 !=0)
    stop("wrong number of params")
  index.sin <- seq(1,lengthParams,by=2)
  # function to compute jacobian of the transformation sinCos2AmplitudeShift()
  one <- function(i=1,params){
    coef.sin <- params[i]
    coef.cos <- params[i+1]

    dAmplitude.dcoef.sin <- coef.sin/sqrt(coef.cos^2+coef.sin^2)
    dAmplitude.dcoef.cos <- coef.cos/sqrt(coef.cos^2+coef.sin^2)

    dShift.dcoef.sin <- - coef.cos/(coef.cos^2+coef.sin^2)
    dShift.dcoef.cos <- coef.sin/(coef.cos^2+coef.sin^2)
    return(c(dAmplitude.dcoef.sin,dShift.dcoef.sin,dAmplitude.dcoef.cos,dShift.dcoef.cos))
  }
  jacobi<-sapply(index.sin,one,params=params)
  res <- matrix(0,nrow=lengthParams,ncol=lengthParams)
  j<-0
  for (i in index.sin){
    j<-j+1
    res[i:(i+1),i:(i+1)] <- jacobi[,j]
  }
  return(res)
}



###################################################
### chunk number 12: 
###################################################
## additional (undocumented) functions needed for algo.hhh

######################################################################
# Function to unpack params and ensure that autoregressive parameters 
# lambda and phi are positive
# and proportion parameter is 0 < pi < 1
#
# theta -  (alpha_i, lambda, phi, prop, beta, gamma_i, delta_i, psi)
# designRes -  result of a call to make.design
######################################################################
unpackParams <- function(theta, designRes){
  
  dimIntercept <- designRes$dimTheta$intercept
  dimLambda <- designRes$dimTheta$lambda
  indexLambda <- dimIntercept+dimLambda
  dimPhi <- designRes$dimTheta$phi
  indexPhi <- indexLambda +dimPhi
  dimProportion <-  designRes$dimTheta$proportion
  indexProportion <- indexPhi+dimProportion
  dimTrend <- designRes$dimTheta$trend
  indexTrend <- indexProportion+dimTrend
  dimSeason <- designRes$dimTheta$season
  indexSeason <- indexTrend +dimSeason
  dimNegbin <- designRes$dimTheta$negbin

  # params set to NULL if not specified
  # intercept always
  alpha <- theta[1:dimIntercept]
  
  if(dimLambda >0)
     lambda <- exp(theta[(dimIntercept+1):indexLambda])
  else lambda <- NULL
  
  if(dimPhi >0)
    phi <- exp(theta[(indexLambda+1):(indexPhi)])
  else phi <- NULL
  
  if(dimProportion >0){
    prop <- theta[(indexPhi+1):indexProportion]
    # ensure that proportion is 0<pi<1
    prop <- exp(prop)/(1+exp(prop))
  } else prop <- NULL
  
  if(dimTrend >0)
    beta <- theta[(indexProportion+1):indexTrend]
  else beta <- NULL
  
  if(dimSeason >0)
    gamma <- theta[(indexTrend+1):indexSeason]
  else  gamma <- NULL
  
  if(dimNegbin >0)
    psi <- exp(theta[(indexSeason+1):(indexSeason+dimNegbin)])
  else psi <- NULL
  
  return(list(alpha=alpha,lambda=lambda, phi=phi,pi=prop,beta=beta, gamma=gamma, psi=psi))
}


#############################################
# function to compute gradient of loglikelihood
# -> used in optim
################################################
gradient <- function(theta,designRes){
  
  if(any(is.na(theta) | !is.finite(theta))) 
    return(rep(NA,length(theta)))
  
  Y<-designRes$Y
  Ym1 <-designRes$Ym1
  control <- designRes$control
  
  mean <- meanResponse(theta=theta, designRes=designRes)
  params <- unpackParams(theta,designRes)
  nOfNeighbours <- designRes$nOfNeighbours
  nhood <- designRes$disProgObj$neighbourhood
  nareas <- ncol(Y)
    
  endemic <- mean$endemic
  meanTotal <- mean$mean
  

  ## helper function for derivatives:
  # negbin model or poisson model
  if(control$negbin!="none"){
    psi <- matrix(params$psi,ncol=nareas,nrow=nrow(Y),byrow=TRUE)
    psiPlusMu <- psi + meanTotal
    
    # helper function for derivatives: negbin
    derivHHH <- function(dmu){
#      if(any(dim(dmu)!=dim(Y))) 
#        cat("warning: dimensions wrong \n")
      (-psi/psiPlusMu +Y/meanTotal -Y/psiPlusMu)*dmu
    }
      
  } else {
    # helper function for derivatives: poisson
    derivHHH <- function(dmu){
#      if(any(dim(dmu)!=dim(dmu))) 
#        cat("warning: dimensions wrong \n")
      Y *(dmu/meanTotal) - dmu
    }

  }

  ###########################################
  ## epidemic part
  ##########################################
  
  # model with lambda and phi
  if(designRes$dimTheta$proportion == 0){
    # gradient for lambda
    if(designRes$dimTheta$lambda >0){
      lambda <- params$lambda
      
      if(length(control$lambda)>1){
        # create vector lambda with elements 0 if control$lambda=FALSE
        lambda <- rep(0,nareas)
        lambda[!is.na(designRes$control$lambda)] <- params$lambda
      }
      
      lambda <- matrix(lambda,ncol=nareas,nrow=nrow(Y),byrow=TRUE)
      dLambda <- derivHHH(lambda*designRes$Ym1)

      # multiple lambda_i's or single lambda ?
      if(length(control$lambda) > 1)
        grLambda <- colSums(dLambda)[!is.na(designRes$control$lambda)]
      else grLambda <- sum(dLambda)

      if(any(is.na(grLambda))){
        warning("derivatives for lambda not computable\n")
        return(rep(NA,length(theta)))
      }

    } else grLambda <- NULL
    

    # gradient for phi
    if(designRes$dimTheta$phi >0){
      phi <- params$phi
      
      if(length(control$neighbours)>1){
        # create vector phi with elements 0 if control$neighbours=FALSE
        phi <- rep(0,nareas)
        phi[!is.na(designRes$control$neighbours)] <- params$phi
      }

      phi <- matrix(phi,ncol=nareas,nrow=nrow(Y),byrow=TRUE)
      if(any(is.na(phi))) stop("phi contains NA\'s\n")
      dPhi <- derivHHH(phi*designRes$Ym1.neighbours)

      # multiple phi_i's or single phi ?
      if(length(control$neighbours)>1)
        grPhi <- colSums(dPhi)[!is.na(designRes$control$neighbours)]
      else grPhi<- sum(dPhi)
      
      if(any(is.na(grPhi))){
        warning("derivatives for phi not computable\n")
        return(rep(NA,length(theta)))
      }

    } else grPhi <- NULL

    # gradient for proportion pi
    grPi <- NULL
    
  } else {
    ################################################
    ## model with lambda and proportion pi
    ###############################################
    
    ## gradient for lambda
    gradLambda <- function(prop,lambda){
    
      # ensure region id is not included
      diag(nhood) <- 0
      # compute lambda_id* [pi_id*Ym1_id + sum_j~id {(1-pi_id )/|j~id|* Ym1_id}] for unit id
      
      dLambda.id <- function(id){
        # number of Neigbours for unit id,  i.e. |k~id| 
        n<-nOfNeighbours[id]
        lambdaYm1.id <- Ym1[,id]*lambda[id]
        pi.id.j <- rep(0,nareas)
        pi.id.j[id]<- prop[id]
        pi.id.j[nhood[,id]>0] <-(1-prop[id])/n
        
        lambdaYm1pi.id <-lambdaYm1.id*matrix(pi.id.j,ncol=nareas,nrow=nrow(Ym1),byrow=TRUE)
        
        # d/dpi log L(mu_i,t)
        return(rowSums(derivHHH(lambdaYm1pi.id)))
      }
      return(sapply(1:nareas,dLambda.id))
    }

    ## gradient for pi
    gradPi <- function(prop,lambda){
    
      # ensure region id is not included
      diag(nhood) <- 0
      # compute (pi_id-pi_id^2)* [lambda_id*Ym1_id - sum_j~id {lambda_id/|j~id|* Ym1_id}] for unit id
      
      dPi.id <- function(id){
        # number of Neigbours for unit id,  i.e. |k~id| 
        n<-nOfNeighbours[id]
        dPiYm1.id <- Ym1[,id]*(prop[id]-prop[id]^2)
        lambda.id.j <- rep(0,nareas)
        lambda.id.j[id]<- lambda[id]
        lambda.id.j[nhood[,id]>0] <-(-lambda[id])/n
        
        dPiYm1lambda.id <-dPiYm1.id*matrix(lambda.id.j,ncol=nareas,nrow=nrow(Ym1),byrow=TRUE)
        
        # d/dpi log L(mu_i,t)
        return(rowSums(derivHHH(dPiYm1lambda.id)))
      }
      return(sapply(1:nareas,dPi.id))
    }

    # gradient for lambda
    if(designRes$dimTheta$lambda ==0)
     cat("no lambda\n")
   
    lambda <- rep(params$lambda,length=nareas)
    prop <- rep(params$pi, length=nareas)

    dLambda <- gradLambda(prop=prop,lambda=lambda)

    # multiple lambda_i's or single lambda ?
    if(designRes$dimTheta$lambda > 1)
      grLambda <- colSums(dLambda)
    else grLambda <- sum(dLambda)
    
    if(any(is.na(grLambda))){
      warning("derivatives for lambda not computable\n")
      return(rep(NA,length(theta)))
    }

    # gradient for phi
    grPhi <- NULL

    # gradient for proportion pi
    dPi <- gradPi(prop=prop,lambda=lambda)
    
    if(designRes$dimTheta$proportion >1)
      grPi <- colSums(dPi)
    else grPi <- sum(dPi)

    if(any(is.na(grPi))){
      warning("derivatives for pi not computable\n")
      return(rep(NA,length(theta)))
    }

  }
  
  
  ############################################
  ## endemic part
  ############################################
  # gradient for intercepts
  grAlpha <- colSums(derivHHH(endemic))
  if(any(is.na(grAlpha))){
    warning("derivatives for alpha not computable\n")
    return(rep(NA,length(theta)))
  }

  # gradient for trend
  if(designRes$dimTheta$trend >0){
    dTrend <- derivHHH(endemic*designRes$X.trendSeason[,1])
    
    if(designRes$dimTheta$trend >1)
      grTrend <- colSums(dTrend)[designRes$control$linear]
    else grTrend <- sum(dTrend)
    
    if(any(is.na(grTrend))){
      warning("derivatives for trend not computable\n")
      return(rep(NA,length(theta)))
    }
    
  } else grTrend <- NULL
  
  # gradient for season
  grSeason <- NULL  
  if(designRes$dimTheta$season >0){

    ## single or multiple seasonal params
    if(length(control$nseason)==1){
    
      for (i in ((designRes$dimTheta$trend>0) +1):ncol(designRes$X.trendSeason) ){
        grSeason <- c(grSeason, sum(derivHHH(endemic*designRes$X.trendSeason[,i])))
      }  
      if(any(is.na(grSeason))){
        warning("derivatives for seasonal parameters not computable\n")
        return(rep(NA,length(theta)))
      }
      
    } else if(length(control$nseason)==nareas){
      #maximum number of Fourier frequencies S.max=max_i{S_i}
      maxSeason <- 2*max(control$nseason)
      
      grSeason <- matrix(NA,nrow=maxSeason,ncol=ncol(Y))
      
      for (j in ((designRes$dimTheta$trend>0) +1):(maxSeason+(designRes$dimTheta$trend>0) ) ){ 
        # compute derivatives of gamma_{ij}, j= 1, ..., 2*S.max
        grSeason[j-(designRes$dimTheta$trend>0),] <- colSums(derivHHH(endemic*designRes$X.trendSeason[,j]))
        # set gradients for gamma_{ij} to NA if  j > S_i
        grSeason[j-(designRes$dimTheta$trend>0),(j > (2*control$nseason)+(designRes$dimTheta$trend>0))] <- NA
      }  
      # gradient now is in order sin(omega_1)_A, sin(omega_1)_B, sin(omega_1)_C, ...  
      #                          cos(omega_1)_A, cos(omega_1)_B, cos(omega_1)_C, ...
      #                          sin(omega_2)_A, sin(omega_2)_B, sin(omega_2)_C, ...  
      #                          ...
      # and needs to be in the following order:
      #     sin(omega_1)_A, cos(omega_1)_A, sin(omega_2)_A, ..., cos(omega_S.max)_A
      #     sin(omega_1)_B, cos(omega_1)_B, sin(omega_2)_B, ..., cos(omega_S.max)_B
  
      # remove NA's, i.e. only derivatives for {gamma_{ij}: j <=2*S_i}
      # check if there are any NaN's
      if(any(is.nan(grSeason))){
        warning("derivatives for seasonal parameters not computable\n")
        return(rep(NA,length(theta)))
      }
      grSeason <- grSeason[!is.na(grSeason)]
 
    } # end multiple params

  } # end gradient season

  # gradient for psi
  if(designRes$dimTheta$negbin>0){
    dPsi <- psi*(digamma(Y+psi)-digamma(psi) +log(psi)+1 - log(psiPlusMu) -psi/psiPlusMu -Y/psiPlusMu)
    
    # multiple psi_i's or single psi?
    if(designRes$dimTheta$negbin >1)
      grPsi <- colSums(dPsi)
    else grPsi <- sum(dPsi)
    
    if(any(is.na(grPsi))){
      warning("derivatives for psi not computable\n")
      return(rep(NA,length(theta)))
    }
    
  } else grPsi <- NULL
  

  res <- c(grAlpha,grLambda,grPhi,grPi,grTrend,grSeason,grPsi)

  return(res)

}

################################
# Calculates the weighted sum of counts of adjacent areas
# weights are specified in neighbourhood-matrix of the disProgObj
# (experimental atm)
# 
# \nu_i,t = \lambda_y_i,t-1 + \phi*\sum_(j~i) [w_ji*y_j,t-1]
#
# disProgObj$neighbourhood can either be a matrix with weights w_ji (in columns)
# or an array (for time varying weights)
#
# if the neighbourhood-matrix has elements 1 if i~j and 0 otherwise
# weightedSumNeighbours() = sumNeighbours()
###########################################
weightedSumNeighbours <- function(disProgObj){

  observed <- disProgObj$observed
  ntime<-nrow(observed)
  narea<-ncol(observed)
  neighbours <- matrix(nrow=ntime,ncol=narea)
  
  nhood <- disProgObj$neighbourhood
  
  #check neighbourhood
  if(any(is.na(nhood)))
    stop("No correct neighbourhood matrix given\n")
  
  ## constant neighbourhood (over time)?
  if(length(dim(nhood))==2){
    # ensure only neighouring areas are summed up
    diag(nhood) <- 0
    nhood <- array(nhood,c(narea,narea,ntime))
  
  } else if(length(dim(nhood))==3){
    if(any(dim(nhood)[1:2]!= narea) | dim(nhood)[3] != ntime) 
      stop("neighbourhood info incorrect\n")
  }
    
  # number of neighbours
  nOfNeighbours <-colSums(nhood[,,1]>0)
    
  for(i in 1:ncol(observed)){
    #weights <- matrix(as.numeric(nhood[,i]),nrow=nrow,ncol=ncol,byrow=TRUE)
    weights <- t(nhood[,i,])
    neighbours[,i] <- rowSums(observed*weights)
  }
  
  return(list(neighbours=neighbours, nOfNeighbours=nOfNeighbours))
}


#################################################
# params psi, lambda and phi are on log-scale
# -> transformation of estimates, standard errors and fisher (using delta rule)
# labels for results
#
# g(theta) = (exp(lambda), exp(phi), beta, gamma, delta, exp(psi), alpha)
# D is the Jacobian of g
# D = diag(exp(lambda), exp(phi), 1, 1, 1, exp(psi), 1)
#########################################
jacobian <- function(thetahat, designRes){  
  
  dimtheta <- designRes$dimTheta$dim
  nareas <- ncol(designRes$disProgObj$observed)
  
  thetaNames <- NULL
  D <-diag(1,ncol=dimtheta,nrow=dimtheta)
  
  dimLambda <- designRes$dimTheta$lambda
  dimPhi <- designRes$dimTheta$phi
  dimPi <- designRes$dimTheta$proportion
  dimTrend <- designRes$dimTheta$trend
  dimPsi <- designRes$dimTheta$negbin
  dimSeason <-designRes$dimTheta$season
  nseason <- designRes$control$nseason
  
  alpha <- colnames(designRes$disProgObj$observed)
  if(is.null(alpha)) alpha <- paste("obs",1:nareas, sep="")
  thetaNames <- c(thetaNames, alpha)

  if(dimLambda >0){
    if(length(designRes$control$lambda)==1)
      lambda <- "lambda"
    else {
      lambda <- paste("lambda", alpha, sep="_")[!is.na(designRes$control$lambda)]
    }
    thetaNames <- c(thetaNames, lambda)
    
    index <-(nareas+1):(nareas+dimLambda)
    thetahat[index] <- exp(thetahat[index])
    diag(D)[index] <- thetahat[index]
  } 
  
  if(dimPhi >0){
    if(length(designRes$control$neighbours)==1)
      phi <- "phi"
    else {
      phi <- paste("phi", alpha, sep="_")[!is.na(designRes$control$neighbours)]
    }
    thetaNames <- c(thetaNames, phi)
    
    index <- (nareas+dimLambda+1):(nareas+dimLambda+dimPhi)
    thetahat[index] <- exp(thetahat[index])
    diag(D)[index] <- thetahat[index]
  }

  if(dimPi>0){
    prop <- switch(designRes$control$proportion,
                    "single"="pi",
                    "multiple"=paste("pi", alpha, sep="_"))

    thetaNames <- c(thetaNames, prop)
    index <- (nareas+dimLambda+dimPhi+1):(nareas+dimLambda+dimPhi+dimPi)
    exp.pi <- exp(thetahat[index])
    diag(D)[index] <- exp.pi/((1+exp.pi)^2)
    thetahat[index] <- exp.pi/(1+exp.pi)
  }
  
  if(dimTrend >0){
    beta <- colnames(designRes$X.trendSeason)[1]
    if(length(designRes$control$linear)>1)
    beta <- paste(beta,alpha,sep="_")[designRes$control$linear]
                   
    thetaNames <- c(thetaNames, beta)
  }
  
  if(dimSeason > 0){
    maxSeason <- 2*max(nseason)
    sinCos <- rep(colnames(designRes$X.trendSeason)[(1+ (dimTrend>0) ):((dimTrend>0) +maxSeason)], length=maxSeason)
    
    if(length(nseason)==1){
      gammaDelta <- sinCos
    } else if(length(nseason==nareas)){
      gammaDelta <- matrix(NA,ncol=nareas,nrow=maxSeason)
      for(i in 1:nareas){
        gammaDelta[0:(2*nseason[i]),i] <- paste(sinCos,alpha[i],sep="_")[0:(2*nseason[i])]
      }
      gammaDelta <- gammaDelta[!is.na(gammaDelta)]
    } 
    thetaNames <- c(thetaNames, gammaDelta )
  }
  if(dimPsi >0){
    psi <- switch(designRes$control$negbin,
                  "single"="psi",
                  "multiple"=paste("psi",alpha,sep="_"))
                  
    thetaNames <- c(thetaNames, psi)
    index <- (dimtheta-dimPsi+1):dimtheta
    thetahat[index] <- exp(thetahat[index])
    diag(D)[index] <-  thetahat[index]
  }

  dimnames(D) <- list(thetaNames,thetaNames)
  names(thetahat) <- thetaNames
  return(list(D=D,theta=thetahat))
}

# theta.epidemic = c(lambda,phi)
# Note: lambda and phi are  on log-scale
getLambda <- function(theta.epidemic, designRes, t.weights=1){
  
  # check dimension of theta.epidemic
  dimLambda <- designRes$dimTheta$lambda
  dimPhi <- designRes$dimTheta$phi
  if(designRes$dimTheta$proportion>0)
    stop("proportions currently not supported\n")
  
  if(length(theta.epidemic)!= (dimLambda+dimPhi))
    stop("vector with parameters must be of length ", dimLambda+dimPhi,"\n")
    
  # is there an autoregression?
  if(sum(!is.na(designRes$control$lambda))==0 & sum(!is.na(designRes$control$neighbours)) ==0)
    return(NULL)

  if(dimLambda>0){
    coef.lambda <- exp(theta.epidemic[1:dimLambda] )
  } else coef.lambda <- 0
  if(dimPhi>0){
    coef.phi <- exp(theta.epidemic[(dimLambda+1):length(theta.epidemic)] )
  } else coef.phi <- 0

  #univariate?
  if(ncol(designRes$disProgObj$observed)==1){
    if(sum(!is.na(designRes$control$lambda))==1)
    return(coef.lambda)
    else return(NULL)
  }
  
   nareas <- ncol(designRes$Y) #ncol(nhood)
  
  if(designRes$control$proportion=="none"){
    # no lambda
    if(sum(!is.na(designRes$control$lambda))==0){
      lambda <- rep(0,nareas)
    # single lambda for all units
    } else if(sum(!is.na(designRes$control$lambda))==1 & length(designRes$control$lambda)==1){
      lambda <- rep(coef.lambda,nareas)
    # multiple lambda
    } else{
      lambda <- rep(0, nareas)
      lambda[designRes$control$lambda] <- coef.lambda
    }
    Lambda <- diag(lambda,nareas)
    
    if(dimPhi>0){
      # extract neighbourhood, i.e. weight matrix
      nhood <- designRes$disProgObj$neighbourhood
      # time-varying weights w_ji
      if(length(dim(nhood))==3)
        nhood <- nhood[,,t.weights]
      
      # ensure the diagonal is zero
      diag(nhood) <- 0
      nOfNeighbours <- colSums(nhood>0)
    
      # single phi for all units
      if(length(designRes$control$neighbours)==1 & sum(!is.na(designRes$control$neighbours))==1){
        phi <-rep(coef.phi,nareas)
      } else if(length(designRes$control$neighbours)>1 & sum(!is.na(designRes$control$neighbours))>0){
        phi <- rep(0,nareas)
        phi[!is.na(designRes$control$neighbours)] <- coef.phi
      }  
      phi.weights <- matrix(phi,nrow=nareas,ncol=nareas,byrow=FALSE)*nhood
      Lambda[nhood>0] <- phi.weights[nhood>0]
    }

  } else { #todo: check
    return(NULL)
#hoehle 14 Oct 2008 - commented, coz it contains warnings for R CMD check
#    lambdaMatrix <- matrix(lambda,ncol=nareas,nrow=nareas,byrow=TRUE)
#    nOfNeighbours <- rowSums(nhood)
#    piMatrix <- matrix((1-prop)/nOfNeighbours,ncol=nareas,nrow=nareas,byrow=TRUE)
#    piMatrix[nhood==0] <-0
#    diag(piMatrix)<-prop
#    Lambda <- lambdaMatrix*piMatrix
  }
  return(Lambda)
}

## moment estimator of exp(alpha)
## alpha.hat(lambda,phi) = mean(y)' %*% (I - Lambda)
expAlpha.mm <- function(Lambda,Y){
  mean.obs <- colMeans(Y)
  mean.obs %*% (diag(1,length(mean.obs))-Lambda)
}


########
logLik.ah <- function(object,...){
   if(!inherits(object, "ah"))
      stop("expected object to be an object of class ah\n")
   if(!object$convergence)
      stop("algorithm did not converge\n")
   val <- object$loglikelihood
   attr(val, "df") <- length(coef(object))
   attr(val, "nobs") <- object$nObs
   class(val) <- "logLik"
   return(val)
}
logLik.ahg <- function(object, ...){
   logLik.ah(object$best)
}



