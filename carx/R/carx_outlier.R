
#' Detecting outliers in a  \code{carx} model.
#'
#' This is an internal function. It tests for additional outlier one at a time, for each time point, as described in
#'  Wang and Chan (2015), adjusted for multiplicity of testing. If the test result is significant, the function 
#'  augments the location of the most significant outlier to the vector of outlier indices, 
#'  i.e, the \code{outlier.indices} in the object returned by the function.
#' @param object a \code{carx} object
#' @keywords internal
#' @return a possibly updated \code{object} which will have an attribute \code{outlier.indices} 
#' denoting the indices of outliers with the new index of outlier appended, if any outlier is detected.
ot.carx <- function(object)
{
	#message("detecting outliers")
	nSample <- 10000
	threshold <- 0.025/object$nObs
	eps <- stats::rnorm(nSample,0,object$sigma)
	trend <- object$x%*%object$prmtrX
	covEta <- computeCovAR(object$prmtrAR, object$sigma,lag=object$p)
	nObs <- object$nObs
	p <- object$p
	prmtrAR <- object$prmtrAR
	skipIndex <- object$skipIndex
	y <- object$y
	ci <- object$ci
	lcl <- object$lcl
	ucl <- object$ucl
  finiteRows <- object$finiteRows

	if(any(ci[finiteRows]>0))
    y[finiteRows][ci[finiteRows]>0] <- ucl[finiteRows][ci[finiteRows]>0]
  if(any(ci[finiteRows]<0))
	y[finiteRows][ci[finiteRows]<0] <- lcl[finiteRows][ci[finiteRows]<0]

	pValues <- numeric(nObs)
	pValues[skipIndex] <- 1

	for(idx in seq(1,nObs)[-skipIndex])
	{
		#message(sprintf("checking %i",idx))
		wkm <- y[(idx-1):(idx-p)]
		tmpCensorIndicator <- ci[(idx-1):(idx-p)]
		nCensored <- sum(tmpCensorIndicator!=0)
    censored <- tmpCensorIndicator[tmpCensorIndicator!=0]
		if(nCensored)   #at least one is censored
		{
			if( nCensored < p ) #not all are censored
			{
				conditionalIndex <- which(tmpCensorIndicator==0) #indices of known values
				tmpY <- y[(idx-1):(idx-p)][conditionalIndex] #known values
				tmpM <- trend[(idx-1):(idx-p)]
				cdist <- conditionalDistMvnorm(tmpY, conditionalIndex,tmpM,covEta)
				tmpMean <- cdist$'mean'
				tmpVar <- cdist$'var'
			}else{
				tmpMean <- trend[(idx-1):(idx-p)]
				tmpVar <- covEta
			}
      tmpLower <- rep(-Inf,length = nCensored)
      tmpUpper <- rep(Inf,length = nCensored)
      tmpLower[censored>0] <- object$ucl[(idx-1):(idx-p)][tmpCensorIndicator>0]
      tmpUpper[censored<0] <- object$lcl[(idx-1):(idx-p)][tmpCensorIndicator<0]
			smpl <- tmvtnorm::rtmvnorm(nSample,tmpMean,tmpVar,lower=tmpLower,upper=tmpUpper,algorithm="gibbs")
			smpl <- as.matrix(smpl)
			ySmpl <- numeric(nSample)
			for(i in 1:nSample){
				wkm[tmpCensorIndicator!=0] <- smpl[i,]
				ySmpl[i] <- trend[idx] + (wkm - trend[(idx-1):(idx-p)])%*%prmtrAR + eps[i]
			}
			pU <- sum(ySmpl > y[idx])/nSample
			pL <- sum(ySmpl < y[idx])/nSample
		}
		else{
			r <- y[idx]-trend[idx] - (wkm-trend[(idx-1):(idx-p)])%*%prmtrAR
			r <- r/object$sigma
			pU <- stats::pnorm(r,lower.tail=FALSE)
			pL <- stats::pnorm(r,lower.tail=TRUE)
		}
		pValues[idx] <- min(pU,pL)
	}
	minP <- min(pValues)
	if( minP <= threshold )
  {
		i <- which(pValues == minP)
    if(is.null(object$outlier.indices))
       object$outlier.indices <- i
    else
      if(!(i %in% object$outlier.indices))
       object$outlier.indices <- c(object$outlier.indices, i)
  }
  object
}

#' S3 method to detect outlier of a \code{carx} object
#'
#' Detect all outliers of a \code{carx} object.
#' @inheritParams outlier.carx
#' @return an updated \code{carx} object. If any outlier is detected, its index will be stored in the \code{outlier.indices} attribute of the return object, and prefix for variable name is stored in the \code{outlier.prefix} attribute. Note that if the original object is fitted through a formula interface, the formula will also be updated.
#' @export
#' @seealso \code{\link{outlier.carx}}.
#'
outlier <- function(object,outlier.prefix="OI_",seed=NULL) UseMethod("outlier")

#' Detect all outliers of a \code{carx} object
#'
#' Detect all outliers of a \code{carx} object and update the model if any outlier is detected.
#' It tests for the presence of outliers one at a time, for each time point, adjusted for multiplicity of testing, as described in Wang and Chan (2015).
#' @param object a \code{carx} object.
#' @param outlier.prefix the prefix used to construct variable name for indicator variables representing the detected outliers, default = "OI_".
#' @param seed the seed for randon number generator, default=\code{NULL}.
#' @return an updated \code{carx} object. If any outlier is detected, its index will be stored in the \code{outlier.indices} attribute of the return object, and prefix for variable name is stored in the \code{outlier.prefix} attribute. Note that if the original object is fitted through a formula interface, the formula will also be updated.
#' @references
#' Wang C, Chan KS (2015). "Quasi-likelihood estimation of a censored autoregressive model with exogenous variables." Submitted.
#
#' @export
#' @examples
#' sigma = 0.6
#' nObs = 100
#' dat = carxSimCenTS(nObs=nObs,sigma=sigma,ucl=Inf)
#' dat$y[as.integer(nObs/2)] = dat$y[as.integer(nObs/2)] + 4*sd(dat$y)
#' mdl <- carx(y~X1+X2-1,data=dat, p=2, CI.compute = FALSE)
#' oc = outlier(mdl)
#' #note the outlier indices in the output:
#' print(oc)
#' #note the updated formula:
#' print(formula(oc))

outlier.carx <- function(object,outlier.prefix="OI_",seed=NULL)
{
  message("Detecting outliers.")
  if(!is.null(seed)) set.seed(seed)
  newObj = object
  newFormula = NULL
  tryCatch({newFormula=formula(object)},error=function(e){newFormula=NULL})

  oiVec = NULL
  ot = -1
  while(TRUE)
  {
    preOt <- newObj$outlier.indices
    newObj = ot.carx(newObj)
    newOt <- newObj$outlier.indices
    if(is.null(newOt))
      break
    if(!is.null(preOt) & length(preOt) == length(newOt))
      break
    oi = numeric(newObj$nObs)
    idx <- newOt[length(newOt)]
    oi[idx] = 1
    newVar = paste0(outlier.prefix,idx)
    if(is.null(newFormula))
    {
      newx=data.frame(newObj$x)
      newx[,newVar] <- oi
      #newObj = update(newObj,x=newx) #this doesn't work.
      newObj=carx(y=newObj$y,
                  x=newx,
                  ci=newObj$ci,
                  lcl=newObj$lcl,
                  ucl=newObj$ucl,
                  p=newObj$p,
                  prmtrX = c(newObj$prmtrX,0),
                  prmtrAR = newObj$prmtrAR,
                  sigma = newObj$sigma,
                  tol = newObj$tol,
                  max.iter = newObj$max.iter,
                  CI.compute=FALSE
                 )
    }
    else
    {
      #we have formula and data

      #update the data
      if(!is.null(newObj$cenTS))
      {
        newData <- newObj$cenTS
        nms <- names(newData)
        #unfortunately I haven't find a way to assign a column to a xts with arbitrary name, hope that
        newData$newVarByChao <- oi
        names(newData) <- c(nms,newVar)
      }
      #update the formula
      newFormula =stats::update(newFormula, stats::as.formula(paste("~.+",newVar)))
      newObj <- carx(newFormula, data = newData,
                   p=newObj$p,
                   prmtrX = c(newObj$prmtrX,0),
                   prmtrAR = newObj$prmtrAR,
                   sigma = newObj$sigma,
                   tol = newObj$tol,
                   max.iter = newObj$max.iter,
                   CI.compute=FALSE
                   )
    }
    newObj$outlier.indices <- newOt
    newObj$outlier.prefix <- outlier.prefix
  }
  message("Detecting outliers. Done.")
  newObj
}
