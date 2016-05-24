
#' Predict a regression model with AR errors
#'
#' @param y response data used to fit the model.
#' @param x covariate matrix  used to fit the model.
#' @param prmtrX regression coefficients for covariates.
#' @param prmtrAR AR coefficients for the regression errors.
#' @param sigma innovation standard deviation.
#' @keywords internal
#' @inheritParams predict.carx

predictARX <- function(y,x,prmtrX,prmtrAR,sigma,n.ahead,newxreg,CI.level=0.95)
{
  nObs <- length(y)
  p <- length(prmtrAR)
  stopifnot( nObs >= p )
  iStart <- nObs-p+1 #use only latest p observations
  nStart <- p
  x <- as.matrix(x)
  newxreg <- as.matrix(newxreg)
	newxreg <- rbind(x[iStart:nObs,,drop=FALSE],newxreg)
	eta <- y[iStart:nObs] - x[iStart:nObs,,drop=FALSE]%*%prmtrX
	eta <- c(eta, rep(0,n.ahead))
  yPred <- c(y[iStart:nObs], rep(0,n.ahead))
  coefs <- matrix(rep(0,n.ahead*n.ahead),nrow=n.ahead,ncol=n.ahead)
  predSE <- rep(0,n.ahead)
  for(i in 1:n.ahead)
  {
    eta[nStart+i] <- sum(eta[(nStart+i-1):(nStart+i-p)]*prmtrAR)
    yPred[nStart+i] <- newxreg[nStart+i,]%*%prmtrX + eta[nStart+i]
    coefs[i,i] <- 1
    if(i>1)
    {
      if( min(i-1,n.ahead-1,p) > 0)
      {
        for( j in 1:min(i-1,n.ahead-1,p) )
          coefs[i,] <- coefs[i,] + prmtrAR[j]*coefs[i-j,]
      }
    }
    predSE[i] <- sigma*sqrt(sum(coefs[i,]^2))
  }
  yPred <- yPred[-(1:p)]
  probs <- c((1-CI.level)/2,(1+CI.level)/2)
  q <- stats::qnorm(probs)
  qntl <- matrix(nrow=n.ahead,ncol=2)
  qntl[,1] <- yPred + predSE*q[1]
  qntl[,2] <- yPred + predSE*q[2]
  list("fit"=as.vector(yPred),"se.fit"=predSE,"ci"=qntl)
}


#' Prediction with a fitted \code{carx} object
#'
#' Predict the future values of an fitted \code{carx} object. If the model has non-null covariate \code{x} other than the constant mean, the new observations in \code{x} must be supplied via \code{newxreg}. The model prediction is done in a similar way as in \code{fitted.carx} by identifying the latest \eqn{p} consecutive observed responses in the data used to estimate model, then compute the mean of the conditional distribution of the future values given the information since the latest p consecutive observed values and the supplied new covariate values. For more details, see Wang and Chan (2015).
#' @references
#' Wang C, Chan KS (2015). "Quasi-likelihood estimation of a censored autoregressive model with exogenous variables." Submitted.
#'
#'
#' @param object a fitted \code{carx} object.
#' @param newxreg the new observations for the coverates \code{x}, default={NULL}.
#' If there is no covariate, the value can be assigned to be \code{NULL}.
#' Otherwise, a matrix of new observations is required to compute predictions. Default=\code{NULL}.
#' @param n.ahead the number of steps ahead to predict, default = 1.
#' @param CI.level the CI.level to construct the confidence interval, default = 0.95.
#' @param nRep the number of replications to be performed in the bootstrap for prediction confidence intervals when censoring exists in the last \code{p}
#' observations, default = 1000.
#' @param na.action how should \code{model.frame} deal with \code{NA} values.
#' @param useSimulatedResidual if True, the innovation sequence will be sampled from the simulated residuals, otherwise, they will be sampled from normal distribution with mean 0 and standard deviation \code{object$sigma}. Set it \code{TRUE} if the normal distribution may be too restrictive. Default = FALSE.
#' @param ... not used.
#' @return A list consisting of \code{fit}, \code{se.fit}, and \code{ci} representing the predictions,
#' standard errors of predictions, and confidence intervals respectively.
#' @export
#' @examples
#'  #This is the function to run a simulation study about the empirical coverage rate of 
#'  #the predictive confidence intervals,
#'  runSimPredCR <- function(nRep=200,nObs=200,n.ahead=10,
#'                           saveRslt=FALSE,saveDir='./testPredictCR',
#'                           seed=NULL)
#'  {
#'    if(!is.null(seed))
#'      set.seed(seed)
#'    crMat = matrix(nrow=n.ahead,ncol=nRep)
#'    if(saveRslt)
#'    {
#'      dir.create(saveDir,showWarnings=FALSE,recursive=TRUE)
#'    }
#'
#'    replication = list()
#'    simSingle <- function(iRep)
#'    {
#'      message(sprintf("iRep: %04d",iRep))
#'      sdata = carxSim(nObs=nObs+n.ahead)
#'      trainingData = sdata[1:nObs,]
#'      testData = sdata[-(1:nObs),]
#'      mdl = carx(y~X1+X2-1,data=trainingData,p=2)
#'      newxreg = testData[,c('X1','X2')]
#'      predVal = predict(mdl,newxreg=newxreg,n.ahead=n.ahead)
#'      crInd = (predVal$ci[,1] <= testData$y) & (predVal$ci[,2] >= testData$y)
#'      crMat[,iRep] = crInd
#'      list(trainingData=trainingData,
#'           testData=testData,
#'           fitted=mdl,
#'           predVal = predVal,
#'           crInd= crInd)
#'    }
#'    replication = lapply(1:nRep,simSingle)
#'    crMat = sapply(replication,FUN=function(x){x$crInd})
#'    print(crMat)
#'    crPred = apply(crMat,1,mean)
#'    message("empirical coverage rate:")
#'    print(crPred)
#'    if(saveRslt)
#'    { 
#'      save(replication,file=paste0(saveDir,'/replication.RData'))
#'      save(crMat,file=paste0(saveDir,'/crMat.RData'))
#'      save(crPred,file=paste0(saveDir,'/crPred.RData'))
#'    }
#'
#'    list(replication=replication,crMat=crMat,meanCR=crPred)
#'  }

#'  #note that nRep=2 is for illustration only, for more stable result, use nRep>=500.
#'  simPredCR = runSimPredCR(nRep=2,nObs=100)


#'  \dontrun{
#'    # for more stable simulation result, run with nRep = 500.
#'    simPredCR = runSimPredCR(nRep=500,nObs=100)
#'    message("Empirical coverage rate:")
#'    print(simPredCR$meanCR)
#'  }

predict.carx <- function(object,newxreg=NULL,n.ahead=1,CI.level=0.95,nRep=1000,na.action=NULL,useSimulatedResidual=FALSE,...)
{

  if( n.ahead < 1)
	  stop("ERROR: n.ahead must be greater than or equal to 1.")

  if(object$xIsOne) #no x, add intercept
  {
    newxreg <- as.matrix(rep(1,n.ahead))
  }
  else
  {
    if(is.null(newxreg))
      stop("ERROR: newxreg supplied is NULL, but the x data in model is not ones.")

    #browser()
    if(is.vector(newxreg))
      newxreg <- as.matrix(newxreg)
    if(dim(newxreg)[1] != n.ahead)
        stop("New data doesn't have the same row of data as n.ahead.")

    #if formula is in the model, the user is more likely to supply the raw data,
    tryCatch({fml=formula(object)},error=function(e){fml=NULL})
    if(is.null(fml))
    {
      #no formula, the new data must have the same dimension as x
      #if(is.vector(newxreg))
      if(dim(newxreg)[2] != dim(object$x)[2])
        stop("New data doesn't have the same number of variables as object$x.")
    }
    else
    {
      message("I am trying to combine the formula and your supplied data.")
      cfml <- fml
      newxreg <- as.data.frame(newxreg)
      if(!is.null(object$outlier.indices))
      {
        #add outlier variables to the data
        for(idx in object$outlier.indices)
        {
          newVar = paste0(object$outlier.prefix,idx)
          nms <- names(newxreg)
          newxreg$newVarByChao <- numeric(n.ahead)
          names(newxreg) <- c(nms,newVar)
        }
      }
      cfml <- nlme::getCovariateFormula(cfml)
      mf <- stats::model.frame(formula=cfml,data=as.data.frame(newxreg),na.action=NULL)
      newxreg <- stats::model.matrix(attr(mf,"terms"),data=mf)
    }
  }

  p <- object$p
  nObs <- object$nObs
  probs <- c((1-CI.level)/2,(1+CI.level)/2)
  qntl <- matrix(nrow=n.ahead,ncol=2)

  if(dim(newxreg)[1] != n.ahead)
    stop("ERROR: number of rows in x doesn't equal to n.ahead.")

	#find the beginning index of p consecutive observations
	iStart <- 1
	for(i in (nObs-p+1):1)
	{
		if(all(object$ci[i:(i+p-1)]==0))
		{
			iStart <- i
			break
		}
	}

  stopifnot(all(is.finite(object$y[iStart:nObs]))) #this means there exists na values from iStart which we cannot handle right now.

	nStart <- nObs - iStart + 1

	if(nStart == p & !useSimulatedResidual)
	{
	  message("Latest p observations are not censored")
    predictARX(object$y,object$x,object$prmtrX,
               object$prmtrAR,object$sigma,
               n.ahead,newxreg,CI.level=CI.level)
  }
	else
	{
	  message("Latest p observations are censored or we'll use simulated residuals")
	  eta <- object$y[iStart:nObs] - object$x[iStart:nObs,]%*%object$prmtrX
  	eta <- c(eta, rep(0,n.ahead))
    yPred <- c(object$y[iStart:nObs], rep(0,n.ahead))
    predSE <- rep(0,n.ahead)
	  newxreg <- rbind(object$x[iStart:nObs,],newxreg)
		tmpCensorIndicator <- object$ci[nObs:iStart] #reverse order
		nCensored <- sum(tmpCensorIndicator!=0)
    covEta <- computeCovAR(object$prmtrAR, object$sigma, nStart)
    trend <- as.vector(newxreg[nStart:1,]%*%object$prmtrX)

    if(nCensored > 0)
    {
      if( nCensored < nStart )
      {
        conditionalIndex <- which(tmpCensorIndicator==0)
        tmpY <- object$y[nObs:iStart][conditionalIndex]
        cdist <- conditionalDistMvnorm(tmpY, conditionalIndex,trend,covEta)
        tmpMean <- cdist$'mean'
        tmpVar <- cdist$'var'
      }else
      {
        tmpMean <- trend
        tmpVar <- covEta
      }
      tmpLower <- rep(-Inf,length = nCensored)
      tmpUpper <- rep(Inf,length = nCensored)
      censored <- tmpCensorIndicator[tmpCensorIndicator!=0]
      tmpLower[censored>0] <- object$ucl[nObs:iStart][tmpCensorIndicator>0]
      tmpUpper[censored<0] <- object$lcl[nObs:iStart][tmpCensorIndicator<0]
      yCensored <- tmvtnorm::rtmvnorm(nRep,tmpMean,tmpVar,lower = tmpLower,upper=tmpUpper)#,algorithm="gibbs")
    }

    if(useSimulatedResidual)
    {
      simRes <- residuals(object,type="raw")
      eps <- matrix(sample(simRes,nRep*n.ahead,replace=T),nrow=nRep,ncol=n.ahead)
    }else
      eps <- matrix(stats::rnorm(nRep*n.ahead,0,object$sigma),nrow=nRep,ncol=n.ahead)

    etaFuture <- matrix(nrow=nRep,ncol=nStart+n.ahead)

    for(iRep in 1:nRep)
    {
      etaFuture[iRep,nStart:1] <- eta[nStart:1]
      if(nCensored>0) #censoring exists, use simulated values for censored y's.
        etaFuture[iRep,nStart:1][tmpCensorIndicator!=0] <- yCensored[iRep,] - trend[tmpCensorIndicator!=0]
      for(i in 1:n.ahead)
        etaFuture[iRep,nStart+i] <- etaFuture[iRep,(nStart+i-1):(nStart+i-p)]%*%object$prmtrAR + eps[iRep,i]
    }
    center <- newxreg[-(1:nStart),]%*%object$prmtrX
    for(i in 1:n.ahead)
	    qntl[i,] <- stats::quantile(etaFuture[,nStart+i],probs=probs)
    qntl[,1] <- qntl[,1] + center
    qntl[,2] <- qntl[,2] + center
    yPred <-  colMeans(etaFuture[,-(1:nStart)])
    yPred <-  yPred + center
    predSE <- matrixStats::colSds(etaFuture[,-(1:nStart)])
    list("fit"=as.vector(yPred),"se.fit"=predSE,"ci"=qntl)
  }
}

#debug(predict.carx)
