
#' @export
carx <- function(y,...) UseMethod("carx")


#' The default estimation method for a CARX model
#'
#' Estimate a CARX model, and compute the standard errors and confidence intervals of the parameter 
#' estimates by parametric bootstrap.
#' 
#' @param y a vector of possibly censored responses. 
#' @param x a matrix of covariates, or some object which can be coerced to matrix. Default=\code{NULL}, 
#'  indicating a pure AR model for \code{y}.
#' @param ci a vector of numeric values. If its i-th value is zero (negative, positive), then the 
#'  corresponding element of 
#' \code{y} is uncensored (left, right censored). 
#'  Default = \code{NULL}, indicating no censoring.
#' @param lcl a vector of lower censoring limits, or a number assuming constant limit. 
#'  Default = \code{NULL}, indicating no lower censoring limit.
#' @param ucl a vector of upper censoring limits, or a number assuming constant limit.
#'  Default = \code{NULL}, indicating no upper censoring limit.
#' @param p the AR order  for the regression errors. Default = 1.
#' @param prmtrX the initial values of the regression parameters for \code{x}.
#'  Default = \code{NULL}.
#' @param prmtrAR the initial values of the AR coefficients.
#'  Default = \code{NULL}.
#' @param sigma the initial value of the innovation (noise) standard deviation.
#'  Default = \code{NULL}.
#' @param y.na.action a string indicating how to deal with missing (NA) values in \code{y}. 
#'  If it is set to "skip", cases containing a missing value will be skipped, so that the 
#'  estimating equation of future cases will be conditional on 
#'  the most recent \code{p} complete cases after the skipped case. 
#'  If "as.censored", the \code{y} value will be 
#'  treated as left-censored with lower censoring limit replaced by positive infinity. 
#'  The user may choose to use "skip" if there exist few long gaps in the time series of
#'  response. Use "as.censored" if the missing values in \code{y} are many and scattered in time. 
#'  N.B.: The presence of any missing values in \code{x} will automatically hard-code \code{y.na.action} 
#'  to be "skip".
#' @param addMu logical, indicating whether to include an intercept in case \code{x=NULL}. 
#'  Default = \code{TRUE}.
#' @param tol the tolerance level used in the stopping criterion. Default = 1.0e-4.
#' @param max.iter maximum number of iterations allowed in the optimization. Default = 100.
#' @param CI.compute logical value to indicate whether to compute the confidence intervals for the
#'  parameters. Default = \code{FALSE}, as it can be time-consuming to run the parametric bootstrap.
#' @param CI.level numeric value in (0,1) representing the confidence level. Default = 0.95.
#' @param b number of bootstrap replicates used for computing the boostrap confidence intervals.
#'  Default = 1000.
#' @param b.robust logical, if \code{TRUE}, the innovations are re-sampled from the simulated residuals; 
#'  otherwise they are re-sampled from a centered normal distribution with the estimated standard deviation.
#'  Default = \code{FALSE}.
#' @param b.show.progress logical, if \code{TRUE}, a text bar will be displayed to show the progress of bootstrap when computing the confidence intervals of the parameter estimates. 
#' @param init.method a string selecting a procedure ("biased", or "consistent") 
#'  for deteriming the initial values, in case there 
#'  are no initial values for some parameters, i.e., one of \code{prmtrX}, \code{prmtrAR}, \code{sigma} 
#'  is \code{NULL}. The "biased" method is always available, 
#'  as it uses maximum gaussian likelihood estimator with the data replacing any censored observation by its 
#'  corresponding censoring limit. The "consistent" method is only available for left censored data. 
#'   Note that the "consistent" method may produce initial estimates with non-stationary AR coefficients 
#'   or negative innovation variance, in which case, the function will fall back to the "biased" method 
#'   for constucting the initial values (when this happens, the attribute \code{fallBackToBiased} of the returned object will 
#'   be set to \code{TRUE}.)  Default="biased".
#' @param cenTS an optional argument to pass a \code{cenTS} object which contain the data used, 
#'  when \code{carx.formula} is invoked. Default = \code{NULL}.
#' @param verbose logical value indicates whether to print intermediate results for monitoring the 
#'  the progress of the iterative estimation procedure. Default = FALSE.
#' @param ... not used.
#' @return a CARX object of the estimated model.
#' @export
#' @examples
#' dat = carxSim(nObs=100,seed=0)
#' mdl <- carx(y=dat$y, x=dat[,c("X1","X2")], ci=dat$ci, lcl=dat$lcl, ucl=dat$ucl, p=2)
#' #or simply call
#' mdl <- carx(y~X1+X2-1,data=dat, p=2, CI.compute = FALSE)

carx.default <- function(y,x=NULL,ci=NULL,lcl=NULL,ucl=NULL,
       p=1,prmtrX=NULL,prmtrAR=NULL,sigma=NULL,
       y.na.action=c("skip","as.censored"),
       addMu=TRUE,
			 tol=1e-4,max.iter=500,CI.compute=FALSE,CI.level=0.95,b=1000,b.robust=FALSE,b.show.progress=FALSE,
       init.method = c("biased","consistent"),
			 cenTS=NULL,verbose=FALSE,...)
{
  #message("Enter carx.default.")
	verbose <- verbose || options()$verbose
	nObs <- length(y)
  
  y.na.action=match.arg(y.na.action)
  nonfiniteYAsCensored = FALSE
  if(y.na.action == "as.censored")
    nonfiniteYAsCensored = TRUE

  init.method=match.arg(init.method)

  finiteY <- which(is.finite(y))

	#standardize censoreIndicator
  if(is.null(ci))
    ci <- rep(0,length(y))
	ci[finiteY][ ci[finiteY]>0 ] <- 1; ci[finiteY][ ci[finiteY]<0 ] <- -1

	if(!is.null(prmtrAR) & length(prmtrAR) != p)
	  stop("initial values of parameter vector of AR is supplied, but its length (",length(prmtrAR),") is not equal to the order p (",p,")")


	if(is.null(lcl)) # no lower censoring
	{
		if(any(ci[finiteY]<0))
			stop("Error in data: lcl is null but there exist left-censored data.")
		else
			lcl = rep(-Inf,nObs)
	} else{
		if(length(lcl) == 1)
		{
			if(!is.nan(lcl))
				lcl <- rep(lcl,nObs)
			else
			{
				warning("lcl is NaN, I will set it to be -Inf.")
				lcl <- rep(-Inf,nObs)
			}
		}
		else #length > 1
		{
			if(length(lcl) != nObs)
				stop("Error: The dimesion of lower censor limit doesn't match that of y.")
		}
	}

	if(is.null(ucl)) #no upper censoring
	{
		if(any(ci[finiteY]>0)) #check
			stop("Error in data: ucl is null but there exist right-censored data.")
		else
			ucl = rep(Inf,nObs)
	} else{
		if(length(ucl) == 1)
		{
			if(!is.nan(ucl))
				ucl <- rep(ucl,nObs)
			else
			{
				warning("ucl is NaN, I will set it to be Inf.")
				ucl <- rep(Inf,nObs)
			}
		}
		else
		{
			if(length(ucl) != nObs)
				stop("Error: The dimesion of upper censor limit doesn't match that of y.")
		}
	}

	if(any(lcl[finiteY] >= ucl[finiteY]))
		stop("Error in censor limit: some lower censor limits are bigger than upper censor limits.")


	noCensor  <- all(ci[finiteY] == 0)
  if(nonfiniteYAsCensored)
  {
    idx <- !is.finite(y)
    if(any(idx))
    {
      ci[idx] <- -1
      lcl[idx] <- Inf
      ucl[idx] <- -Inf #no need to assign to both censoring limits
      noCensor  <- FALSE
      #treat non-finite Y as left-censored
    }

  }



  if(!is.null(x))
	{
		if(!is.matrix(x)) # x may be a vector
			x <- as.matrix(x)
		if(dim(x)[1] != nObs){
			stop(" The dimension of x doesn't match that of y.")
			return(NULL)
		}
		nX <- dim(x)[2]
		x <- x
		if(all(x==1))
			xIsOne <- TRUE
		else
			xIsOne <- FALSE
	}else
	{
    x <- as.matrix(rep(1,nObs))
    nX <- 1
    xIsOne <- TRUE
	}
	xValid <- TRUE
	if(xIsOne & !addMu)
	  xValid <- FALSE

  #check for finite rows in data and construct skipIndex
  if(!is.null(x))
    finiteRows <- (apply(x, 1, function(x){all(is.finite(x))}))
  else
    finiteRows <- rep(TRUE,length(y))
  if(!nonfiniteYAsCensored)
    finiteRows <- finiteRows & is.finite(y)

  skipIndex <- rep(0,length(y))
  skipIndex[1:p] <- 1
  for(i in 1:length(y))
  {
    if(!finiteRows[i])
      skipIndex[i:(i+p)] <- 1
  }
  skipIndex <- which(skipIndex > 0)
  if(verbose) message("Skip index is constructed as ", paste(skipIndex,collapse=', '))
	nSkip <- length(skipIndex)
  #print(head(x))

	ret = list(y = y,
		   x = x,
		   xIsOne = xIsOne,
		   ci = ci,
		   lcl = lcl,
		   ucl = ucl,
       CI.compute = CI.compute,
       tol = tol,
       max.iter = max.iter,
       CI.level = CI.level,
       b = b,
		   skipIndex = skipIndex,
       finiteRows = finiteRows,
       cenTS=cenTS,
		   nonfiniteYAsCensored=nonfiniteYAsCensored,
       verbose = verbose
		   )
	#print(ret)

	#parameters
	#generic
	if(is.null(prmtrX) || is.null(prmtrAR) || is.null(sigma) )
	{
		prmtrX <- numeric(nX)
		prmtrAR <- numeric(p)
		sigma <- numeric(1)
		needInit <- TRUE
	}
	else
		needInit <- FALSE

	#special to store estimated values
	prmtrXEstd <- numeric(nX)
	prmtrAREstd <- numeric(p)
	sigmaEstd <- numeric(1)

	trend <- numeric(nObs)
	covEta <- matrix(nrow=p+1, ncol = p+1)
	wkMean <- matrix(nrow=nObs, ncol = p+1)
	wkCov <- array(0, dim=c(nObs,p+1,p+1))
	res <- numeric(nObs)

  getNPrmtr <- function(){
		return (nX + p +  1)
	}

  pVal <- numeric(getNPrmtr())


	# working space
	resetWK <- function()
	{
		trend <<- numeric(nObs)
		covEta <<- matrix(nrow=p+1, ncol = p+1)
		wkMean <<- matrix(nrow=nObs, ncol = p+1)
		wkCov <<- array(0, dim=c(nObs,p+1,p+1))
		res <<- numeric(nObs)
    pVal <<- numeric(getNPrmtr())

		#initialize wkMean
		for(idx in seq(1,nObs)[-skipIndex])
			wkMean[idx,] <<- y[idx:(idx-p)]
	}
	#-------------------------------


	getPrmtr <- function(){
		ret <- c(prmtrX, prmtrAR, sigma)
		return(ret)
	}

	getCensorRate <- function(){
		return(sum(abs(ci[finiteRows]))*1.0/sum(finiteRows))
	}


	setInitPrmtrForBootstrap <- function()
	{ # set current parameter to be the estimated
		prmtrX   <<-   prmtrXEstd
		prmtrAR   <<-	prmtrAREstd
		sigma  <<-	sigmaEstd
	}

	setEstdPrmtr <- function()
	{ # store the estimated parameter
		prmtrXEstd   <<- prmtrX
		prmtrAREstd   <<- prmtrAR
		sigmaEstd  <<- sigma
	}

	getEstdPrmtr <- function(){
	  if(xValid)
	    ret <- c(prmtrXEstd, prmtrAREstd, sigmaEstd)
	  else
	    ret <- c(prmtrAREstd, sigmaEstd)
		return(ret)
	}

	updateTrend <- function(){
      trend <<- x%*%prmtrX
	}

	updateCovEta <- function(){
		covEta <<- computeCovAR(prmtrAR, sigma)
	}

	#calculate the expectations
	eStep <- function()
	{
		#if(verbose) message("E-step begin.")
		# expection step
		# need to 1) calculate the joint distriubtion of AR process
		#         2) update the trend function/observations
		#         3) calculate the conditional expection & covariances of the censored observations
		updateCovEta()
		updateTrend()
		for(idx in seq(1,nObs)[-skipIndex])
		{
		  #if(idx ==92 )
			tmpCensorIndicator <- ci[idx:(idx-p)]
			nCensored <- sum(abs(tmpCensorIndicator))
			if(nCensored>0)
			{# at least one is censored
				if( nCensored < (p+1) )
				{
					conditionalIndex <- which(tmpCensorIndicator==0) #indices for observed data
					tmpY <- y[idx:(idx-p)][conditionalIndex]
					tmpM <- trend[idx:(idx-p)]
					cdist <- conditionalDistMvnorm(tmpY,conditionalIndex,tmpM,covEta)
					tmpMean <- cdist$'mean'
					tmpVar <- cdist$'var'
				}else{
					tmpMean <- trend[idx:(idx-p)]
					tmpVar <- covEta
				}
				tmpLower <- rep(-Inf,length = nCensored)
				tmpUpper <- rep(Inf,length = nCensored)
				censored <- tmpCensorIndicator[tmpCensorIndicator!=0]

				# lower limit is upper censor limit
				tmpLower[censored>0] <- ucl[idx:(idx-p)][tmpCensorIndicator>0]
				# upper limit is lower censor limit
				tmpUpper[censored<0] <- lcl[idx:(idx-p)][tmpCensorIndicator<0]
				#if(abs(det(tmpVar)) < 0.001) warning("covariance matrix is nearly singluar.")
				ret <- tmvtnorm::mtmvnorm(tmpMean,tmpVar,lower = tmpLower,upper=tmpUpper)
				wkMean[idx,tmpCensorIndicator!=0] <<- ret$'tmean'
				wkCov[idx,tmpCensorIndicator!=0,tmpCensorIndicator!=0] <<- ret$'tvar'
			}
			#if(verbose) message("E-step done.")
		}
	}

	eStepNaive <- function()
	{
		if(verbose) message("entering eStepNaive")
		if(!noCensor)
		{
			for(idx in seq(1,nObs)[-skipIndex])
			{
				tmpCensorIndicator <- ci[idx:(idx-p)]
				wkMean[idx,tmpCensorIndicator>0] <<- ucl[idx:(idx-p)][tmpCensorIndicator>0]
				wkMean[idx,tmpCensorIndicator<0] <<- lcl[idx:(idx-p)][tmpCensorIndicator<0]
				wkCov[idx,,] <<- 0
			}
		}
		if(verbose) message("E-step Naive done.")
	}

	updatePrmtrAR <- function(){
		wkMean2 <- wkMean
		for(i in 1:(p+1)){
			wkMean2[(p+1):nObs,i] <- wkMean2[(p+1):nObs,i] - trend[(p+2-i):(nObs+1-i)]
		}
		mat <- matrix(0,p,p)
		vec <- numeric(p)
		for(idx in seq(1,nObs)[-skipIndex]){
			vec <- vec + wkMean2[idx,1]*wkMean2[idx,-1] + wkCov[idx,-1,1]
			mat <- mat + outer(wkMean2[idx,-1],wkMean2[idx,-1]) + wkCov[idx,-1,-1]
		}

		tmp <- solve(mat,vec)
		if(any(abs(polyroot(c(1,-tmp)))<=1.0)) #tmp is not feasible
		{
			grad <- mat%*%prmtrAR - vec
			tmp <- prmtrAR - grad
			i <- 1
			while(any(abs(polyroot(c(1,-tmp)))<=1.0) && i <= 10)
			{
				tmp <- prmtrAR - grad/2^i
				i <- i+1
			}
			tmp <- as.vector(tmp)
		}
		return(tmp)
	}

	updatePrmtrEV <- function()
	{
    #tryCatch( {
      if(p > 1){
        tmpY <- wkMean[,1] - wkMean[,-1]%*%prmtrAR
      }
      else{
        tmpY <- wkMean[,1] - wkMean[,-1]*prmtrAR
      }
      tmpY <- tmpY[-skipIndex]

      tmpX <- matrix(0,nObs-nSkip,nX)
      index <- seq(1,nObs)[-skipIndex]
      for(idx in 1:(nObs-nSkip)){
        if (p >1){
          tmpX[idx,] <- x[(index[idx]),] - prmtrAR %*% x[(index[idx]-1):(index[idx]-p),]
        }
        else{
          tmpX[idx,] <- x[(index[idx]),] - prmtrAR * x[(index[idx]-1):(index[idx]-p),]
        }
      }
      ret <- solve(t(tmpX) %*% tmpX, t(tmpX)%*%tmpY )
      #ret <- MASS::ginv(t(tmpX) %*% tmpX) %*%(t(tmpX)%*%tmpY)

      ret
    #},warning=function(w){
    #  message("WARNING:", w)
    #  prmtrX
    #},error=function(e){
    #  message("ERROR:",e)
    #  prmtrX
    #},finally={
    #}
    #)
	}

	updateSigma <- function(){
		idx <- seq(1,nObs)[-skipIndex]
		vec <- wkMean[idx,1] - trend[idx]
		for(i in 1:p){
			vec <- vec - prmtrAR[i]*(wkMean[idx,i+1] - trend[(idx-i)])
		}
		ret <- sum(vec^2)
		tmpVec <- c(1,-prmtrAR)
		for(i in idx){
			ret <- ret + as.numeric( t(tmpVec)%*%wkCov[i,,]%*%tmpVec )
		}
		return( sqrt(ret/(nObs-nSkip)) )
	}

	mStep <- function(){

    prevPrmtr <- c(prmtrX,prmtrAR,sigma)
		delta <- 0

		newPrmtrAR <- updatePrmtrAR()
		delta <- delta + sum((newPrmtrAR - prmtrAR)^2)
		prmtrAR <<- newPrmtrAR

    if(xValid)
    {
      newPrmtrEV <- updatePrmtrEV()
      delta <- delta + sum((newPrmtrEV - prmtrX)^2)
      prmtrX <<- newPrmtrEV
    }

		newSigma <- updateSigma()
		delta <- delta + (newSigma - sigma)^2
		sigma <<- newSigma

    return(sqrt(delta/sum(prevPrmtr)^2))
	}

	estimatePrmtr <- function(tol,max.iter)
	{
		if(verbose) message("estimating parameters")
		delta <- 1.0
		nIter <- 1
		while(delta > tol && nIter < max.iter)
		{
			eStep()
			delta <- mStep()
			if(verbose) message(sprintf("Iter:%i, relDif: %f, prmtr: %s", nIter, delta, toString(round(getPrmtr(),digits = 4)) ) )
			nIter <- nIter + 1
		}
		if(verbose) message(sprintf("Iter:%i, relDif: %f, prmtr: %s \nParameter estimated\n", nIter-1, delta, toString(round(getPrmtr(),digits=4)) ) )
	}

	initPrmtr <- function(tol,max.iter)
	{
		if(verbose) message("Initializing parameters")
		delta <- 1.0
		nIter <- 1

		prmtrX <<- numeric(nX)
		prmtrAR <<- numeric(p)
		sigma <<- 1

		eStepNaive()

		#tmporaly modify skipIndex
		prevSkipIndex <- skipIndex
		skipIndex <<- sort(unique(c(skipIndex,which(!is.finite(rowMeans(wkMean))))))
		nSkip <<- length(skipIndex)
		while(delta > tol && nIter < max.iter)
		{
			updateTrend() # needed since it is done to update the trend which is done in the eStep for our method
			delta <- mStep()
			if(verbose) message(sprintf("Iter:%i, relDif: %f, prmtr: %s", nIter, delta, toString(round(getPrmtr(),digits=4) ) ))
			nIter <- nIter + 1
		}
		#if(verbose) message(sprintf("Iter:%i, relDif: %f, prmtr: %s \nParameter initialized\n", nIter-1, delta, toString(round(getPrmtr(),digits=4) ) ))
    if(verbose) message("Initialize parameters, done.")
		skipIndex <<- skipIndex
		nSkip <<- length(skipIndex)
	}

  initPrmtr2 <- function()
  {
    if(verbose)
      message("Initialize parameters by consistent estimator.")
    #message("Use least squares with instrumental variables to initialize parameters.")
    if(any(is.finite(ucl)))
       warning("Method cannot be used for data with right censoring.")

    U <- NULL #record t such that y[t:(t-p)] are uncensored
    for(t in (1:nObs)[-skipIndex])
    {
      if(all(ci[t:(t-p)]==0))
        U <- c(U,t)
    }
    if(is.null(U))
      message("No consecutive p+1 uncensored observations in the series.")

    nU <- length(U)
    if(verbose) message("# of p+1 uncensored obs:", nU)
    if(!xIsOne)
    {
      z <- matrix(nrow=nU,ncol=nX*(p+1)+p)
      g <- matrix(nrow=nU,ncol=nX*(p+1)+p+1)
      ghat <- matrix(nrow=nU,ncol=nX*(p+1)+p+1)
      for(i in 1:nU)
      {
        t <- U[i]
        z[i,] <- c(x[t,],y[(t-1):(t-p)],as.vector(t(x[(t-(1:p)),])))
      }

      uY <- y[U]
      uYhat <- z%*%solve(t(z)%*%z,t(z)%*%uY)

      h <- uY^2 - lcl[U]*uY

      for( i in 1:nU)
      {
        g[i,] <- c((uY[i] - lcl[U[i]])*z[i,], 1)
        ghat[i,] <- c((uYhat[i] - lcl[U[i]])*z[i,], 1)
      }
      tmp <- solve(t(ghat)%*%g, t(ghat)%*%h)
      prmtrX <<- tmp[1:nX]
      prmtrAR <<- tmp[(nX+1):(nX+p)]
      if(tmp[length(tmp)]>=0)
        sigma <<- sqrt(tmp[length(tmp)])
      else
        sigma <<- NaN
    }
    else
    {
      stop("Not Implemented! for x=1")
    }
    if(verbose) message("Initialize parameters, done.")
  }

	exptdLogLik <- function()
	{
		val <- -(nObs - nSkip)*(1+log(sigmaEstd^2))/2
		val
	}



	bootstrapSample <- function()
	{
    if(b.robust)
    {
      epsPool <- residuals.carx(ret,type="raw")
      eps <- sample(epsPool,nObs)
    }else
      eps <- stats::rnorm(nObs,0, sigma)

		updateTrend()

		eta <- numeric(nObs)
		eta[skipIndex] <- eps[skipIndex]
		for(i in seq(1,nObs)[-skipIndex]){
			eta[i] <- eta[(i-1):(i-p)] %*% prmtrAR + eps[i]
			y[i] <<- trend[i] + eta[i]
		}
		ci <<- rep(0,nObs)
		ci[y<lcl] <<- -1
		ci[y>ucl] <<- 1
		if(verbose) message(paste0("\ncensor rate: ", sum(abs(ci))/nObs))
	}

	bootstrapCI <- function(CI.level,b)
	{

		if(verbose) message(sprintf('Bootstrapping CI'))
		yOriginal <- y #copy y as it will be overwritten


		tmpResult <- matrix(0,b,getNPrmtr())
		if(b.show.progress)pb <- txtProgressBar(1,b,style=3)
		for(i in 1:b){
			#message(sprintf('Bootstrapping CI %i/%i',i,b))
			if(b.show.progress) setTxtProgressBar(pb,i)
			setInitPrmtrForBootstrap()
			bootstrapSample()
			resetWK()
			estimatePrmtr(tol,max.iter)
			tmpResult[i,] <- getPrmtr()
		}
		ret$prmtrEstdBootstrap <<- tmpResult
		if(b.show.progress) close(pb)
		qv <- c((1-CI.level)/2,1-(1-CI.level)/2)
		for(j in 1:(getNPrmtr())){
			CI[j,] <<- stats::quantile(tmpResult[,j],qv)
      pval <- mean(tmpResult[,j]>0)
      pVal[j] <<- 2*min(pval,1-pval)
		}
		covMat <<- stats::cov(tmpResult)
		y <<- yOriginal #set back observed data
	}

	# begin execution
	resetWK()
  fallBackToBiased <- FALSE
	#message("initializing parameters")
	if(needInit || noCensor )  # If no censoring, it is estimating an AR-X
  {
    if(init.method == "biased" | any(ci>0))
      initPrmtr(tol, max.iter)
    else
    {
      if(init.method == "consistent")
        initPrmtr2()
      if(!isStationaryAR(prmtrAR) | is.nan(sigma))
      {
        if(verbose) message("Initial consistent estimator for AR parameter is not stationary, fall back to unbiased estimator")
        initPrmtr(tol, max.iter)
        fallBackToBiased <- TRUE
      }
    }
  }
	prmtrInit <- getPrmtr()
  if(verbose) message("Initial estimator: ",toString(round(prmtrInit,digits=4)))

  #message("skipping true estimation!")
  if(!noCensor)
    estimatePrmtr(tol,max.iter)
	setEstdPrmtr()


	if(xValid)
	  coeff <- c(prmtrXEstd,prmtrAREstd)
	else
	  coeff <- prmtrAREstd
	if(xValid)
	{
	  xnames <- colnames(x)
	  if( is.null(xnames) )
	  {
  	  if(xIsOne)
  	    xnames <- c("mu")
    	else
    		xnames <- paste0('X',1:dim(x)[2])
	  }
	}
	else
	  xnames <- NULL
	names(coeff) <- c(xnames,paste0('AR',1:p))

	ret$coefficients = coeff
	ret$prmtrX = prmtrXEstd
	ret$prmtrAR = prmtrAREstd
	ret$sigma = sigmaEstd
  ret$fallBackToBiased = fallBackToBiased

	ret$censorRate = getCensorRate()
	rnames <- c(names(coeff),"sigma")
	ret$prmtrInit = prmtrInit
	ret$prmtrEstd = getEstdPrmtr()
	names(ret$prmtrInit) <- rnames
	names(ret$prmtrEstd) <- rnames
	ret$nObs = nObs
	ret$logLik = exptdLogLik()
	ret$p = p
	ret$nX = nX
	ret$npar = getNPrmtr()
	ret$aic = 2*(-ret$logLik + ret$npar)
	ret$call = match.call()

	if(CI.compute){
		#rnames <- c(names(coeff),"sigma")
		CI <- matrix(nrow=getNPrmtr(),ncol=2)
		covMat <- matrix(nrow=getNPrmtr(),ncol=getNPrmtr())
		bootstrapCI(CI.level,b)
		rownames(CI) <- rnames
		colnames(CI) <- c(sprintf("%4.2f%%",100*(1-CI.level)/2), sprintf("%4.2f%%",100*(1+CI.level)/2))
		rownames(covMat) <- rnames
		colnames(covMat) <- rnames
		ret$CI <- CI
		ret$vcov <- covMat
    ret$pVal <- pVal
		colnames(ret$prmtrEstdBootstrap) = rnames
	}else{
		ret$CI <- NULL
		ret$vcov <- NULL
    ret$pVal <- NULL
	}
	class(ret) <- "carx"
  #message("Exit carx.default.")
	invisible(ret)
}

#' A formula interface to the \code{carx} method
#'
#' This interface will use the supplied \code{formula} and data provided by \code{data} and other arguments 
#' in \code{...} to invoke the \code{carx.default} method. This is the preferred way of calling the 
#' \code{carx.default} function. The data and censoring information are best passed to the function via
#' a \code{cenTS} object comprising 
#' \code{y}, \code{x}, \code{ci}, \code{lcl}, and \code{ucl}.
#'
#'
#' @seealso \code{\link{cenTS}} on how to construct a \code{cenTS} object.
#' @param formula the formula.
#' @param data the data consisting of the \code{y}, \code{x}, \code{lcl},\code{ci}, and \code{ucl} objects;
#'  can be a \code{list}, \code{data.frame}, or a \code{\link{cenTS}} object.
#' @param ... other arguments supplied to \code{\link{carx.default}}.
#' @export
#' @examples
#' dat = carxSim(nObs=100,seed=0)
#' mdl <- carx(y~X1+X2-1,data=dat, p=2, CI.compute = FALSE)

carx.formula <- function(formula, data=list(),...)
{
  #message("Enter carx.formula")
  vars <- list(...)
  nvars <- names(vars)

  if('cenTS' %in% class(data))
  {
    data2 <- data.frame(zoo::coredata(data))
    vars$cenTS <- data
  }
  else
    data2 <- data

	mf <- stats::model.frame(formula=formula,data=data2,na.action=NULL)
	y <- stats::model.response(mf)
	x <- stats::model.matrix(attr(mf,"terms"),data=mf)

  if( "lcl" %in% names(data2) )
    vars$lcl <- data2$lcl

  if( "ucl" %in% names(data2) )
    vars$ucl <- data2$ucl

  if( "ci" %in% names(data2) )
    vars$ci <- data2$ci

  toPass <- c(list(y=y,x=x),vars)
	#est <- carx.default(y,x,...)
  est <- do.call(carx.default,toPass)
	est$call <- match.call()
  #print(est$call)
	est$formula <- formula
  est$data <- data
  #message("Exit carx.formula")
	est
}

#' The quasi-log-likelihood of a \code{carx} object
#' @param object a fitted \code{carx} object.
#' @param ... not used.
#' @return the quasi-log-likelihood value.
#' @export
#' @examples
#' dat = carxSim(nObs=100,seed=0)
#' mdl <- carx(y~X1+X2-1,data=dat, p=2, CI.compute = FALSE)
#' lk = logLik(mdl)


logLik.carx <- function(object,...)
{
	ret <- object$logLik
	class(ret) <- 'logLik.carx'
	ret
}

#' Compute the  AIC of a fitted \code{carx} object
#'
#' Return the AIC of a fitted \code{carx} object where the maximum log-likelihood is replaced by the 
#' maximum quasi-log-likelihood.
#' @param object a fitted  \code{carx} object.
#' @param ... not used.
#' @param k penalty multiplier for the number of parameters. Default = 2.
#' @return the AIC value = -2*maximum quasi-log-likelihood+k*number of parameters
#' @export
#' @examples
#' dat = carxSim(nObs=100,seed=0)
#' mdl <- carx(y~X1+X2-1,data=dat, p=2, CI.compute = FALSE)
#' ic = AIC(mdl)

AIC.carx <- function(object,...,k=2)
{
	val <- -2*object$logLik + k*object$npar
	#class(val) <- "AIC.carx"
	val
}



#' Print a short description of the fitted model
#' @param x a fitted model object.
#' @param ... not used.
#' @return none.
#' @export
#' @examples
#' dat = carxSim(nObs=100,seed=0)
#' mdl <- carx(y~X1+X2-1,data=dat, p=2, CI.compute = FALSE)
#' print(mdl)


print.carx <- function(x,...)
{
	cat("Call:\n")
	print(x$call)
	cat("\nCoefficients:\n")
	print(x$coefficients)

	cat("\nResidual (innovation) standard deviation:\n")
	print(x$sigma)

	cat("\nCensoring rate:\n")
	print(x$censorRate)
	cat("\nSample size:\n")
	print(x$nObs)

	cat("\nNumber of parameters:\n")
	print(x$npar)

	cat("\nQuasi-log-likelihood:\n")
	print(x$logLik)

	cat("\nAIC:\n")
	print(x$aic)

  if(!is.null(x$outlier.indices))
  {
    cat('\nOutlier Indices:\n')
    print(x$outlier.indices)
  }


	#cat("\nInitial estimates:\n")
	#print(x$prmtrInit)
	#cat("\nFinal estimates:\n")
	#print(x$prmtrEstd)
	if(x$CI.compute){
		cat(paste0("\nConfidence interval:\n"))
		print(x$CI)
		cat("\nVariance-covariance matrix:\n")
		print(x$vcov)
    cat("N.B.: Confidence intervals and variance-covariance matrix\n")
    cat(paste0("are based on ", x$b, " bootstrap samples.\n"))
	}
}

#' Summarize the fitted \code{carx} object
#' @param object a fitted \code{carx} object.
#' @param ... not used.
#' @return a summary.
#' @export
#' @examples
#' dat = carxSim(nObs=100,seed=0)
#' mdl <- carx(y~X1+X2-1,data=dat, p=2, CI.compute = FALSE)
#' summary(mdl)

summary.carx <- function(object,...)
{
	numDig <- function(x,d1,d2){
		y <- x
		y[abs(x) > 0.1^d1] <- round(x[abs(x) > 0.1^d1],d1)
		y[abs(x) < 0.1^d1] <- round(x[abs(x) < 0.1^d1],d2)
		y
	}

	if(is.null(object$CI))
  {
		tab <- cbind(Estimate = object$prmtrEstd)
	}else
  {
		se <- sqrt(diag(object$vcov))
		tVal <- c(coef(object),object$sigma)/se

		tab <- cbind(Estimate = object$prmtrEstd,
			     StdErr =  se[1:(object$npar)],
			     lowerCI = object$CI[1:(object$npar),1],
			     upperCI = object$CI[1:(object$npar),2],
           p.value = object$pVal[1:(object$npar)]
			     )
	}
	res <- list(call=object$call,
		    coefficients=tab,
		    aic = AIC.carx(object))
	class(res) <- "summary.carx"
	res
}

#' Print a summary of an \code{carx} object
#' @param x a summary of an \code{carx} object.
#' @param ... not used.
#' @return none.
#' @export
#' @examples
#' dat = carxSim(nObs=100,seed=0)
#' mdl <- carx(y~X1+X2-1,data=dat, p=2, CI.compute = FALSE)
#' summary(mdl)


print.summary.carx <- function(x,...)
{
	cat("Call:\n")
	print(x$call)
	cat("\n")

	cat("\nCoefficients:\n")
	#print(x$coefficients)
  pval <- 'p.value' %in% colnames(x$coefficients)
  stats::printCoefmat(x$coefficients, P.values = pval, has.Pvalue = pval)
 cat(sprintf("Note: confidence intervals are based on %i bootstraps.  P.values are one-sided.\n", x$b))

	cat("\nAIC:\n")
	print(x$aic)
}
