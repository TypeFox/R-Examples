#' Residuals of a fitted \code{carx} object
#'
#' Computes the residuals of fitted \code{carx} object.
#' When no censoring is present, the ordinary residuals will be computed.
#' Otherwise, the simulated residuals (Gourieroux, Monfort, Renault, and Trognon 1987) of a 
#' fitted \code{carx} object will be computed, as suggested in Wang and Chan (2015).
#' 
#'  The simulated residuals are constructed as follows.
#'  First, impute each unobserved  \eqn{Y_t^*} by a (random) realization from the conditional distribution 
#'  \eqn{D(Y_t^*|\{(Y_s,X_s)\}_{s=1}^t)}, evaluated at the parameter estimate.
#'  Then, refit the model with \eqn{(Y_t^* , X_t)} so obtained, via the method of conditional maximum likelihood; 
#'  the residuals from the latter model are the simulated residuals \eqn{\varepsilon_t}.
#' @references Gourieroux C, Monfort A, Renault E, Trognon A (1987). "Simulated residuals." 
#' Journal of Econometrics, 34(1), 201-252.
#'
#' Wang C, Chan KS (2015). "Quasi-likelihood estimation of a censored autoregressive model with exogenous variables." Submitted.
#' @param object a fitted \code{carx} object.
#' @param type a string indicates which type of residual is to be returned.
#' "raw" returns the (simulated) residuals;
#' "pearson" returns the raw residuals divided by estimated standard error of the residuals.
#' @param seed the seed for the random number generator.
#' @param ... not used.
#' @return the simulated residuals.
#' @export
#' @examples
#' dat = carxSim(nObs=100,seed=0)
#' mdl <- carx(y~X1+X2-1,data=dat, p=2, CI.compute = FALSE)
#' #compute the raw residuals
#' res = residuals(mdl,type="raw")
#' #compute the Pearson residuals
#' res = residuals(mdl,type="pearson")
residuals.carx <- function(object,type=c("raw","pearson"),seed=NULL,...)
{
  type <- match.arg(type)
  if(!is.null(seed)) set.seed(seed)
	nObs <- object$nObs
	p <- object$p
	y <- object$y
  finiteRows <- object$finiteRows

	trend <- object$x%*%object$prmtrX
	eta <- object$y - trend
	for(idx in 1:p)
	{
   if(!finiteRows[idx])
      next
		if(object$ci[idx] > 0)
			y[idx] = object$ucl[idx]
		else if(object$ci[idx] < 0 )
				y[idx] = object$lcl[idx]
	}

	for(idx in (p+1):nObs)
	{
    if(!finiteRows[idx])
      next
		if(object$ci[idx] == 0)
		{
			#message(sprintf("idx %i, not censored",idx))
			next
		}

		if(all(object$ci[(idx-1):(idx-p)]==0))
    {
		  #message(sprintf("idx %i, fast ",idx))
			tmpMean <- trend[idx] + eta[(idx-1):(idx-p)]%*%object$prmtrAR
			if(object$ci[idx] > 0)
				y[idx] <- tmvtnorm::rtmvnorm(1, mean=c(tmpMean), sigma=c(object$sigma),lower=object$ucl[idx],upper=Inf,algorithm="gibbs")
			else
				y[idx] <- tmvtnorm::rtmvnorm(1, mean=c(tmpMean), sigma=c(object$sigma),lower=-Inf, upper = object$lcl[idx],algorithm="gibbs")
			next
		}

		#message(sprintf("idx %i, slow",idx))
		iStart <- 1
		for(i in (idx-p):1)
		{
			if(all(object$ci[i:(i+p-1)]==0))
			{
				iStart <- i
				break
			}
		}

		nStart <- idx - iStart + 1
		#message(sprintf("idx: %i, iStart: %i, nStart: %i",idx,iStart,nStart))
		tmpCensorIndicator <- object$ci[idx:iStart] #reverse order
		nCensored <- sum(tmpCensorIndicator!=0)
		covEta <- computeCovAR(object$prmtrAR, object$sigma, nStart)
		if( nCensored < nStart )
		{
			conditionalIndex <- which(tmpCensorIndicator==0)
			tmpY <- object$y[idx:iStart][conditionalIndex]
			cdist <- conditionalDistMvnorm(tmpY, conditionalIndex, trend[idx:iStart], covEta)
			tmpMean <- cdist$'mean'
			tmpVar <- cdist$'var'
		}else
		{
			tmpMean <- trend[idx:iStart]
			tmpVar <- covEta
		}
		tmpLower <- rep(-Inf,length = nCensored) #( y[idx], censored obs)
		tmpUpper <- rep(Inf,length = nCensored)
		censored <- tmpCensorIndicator[tmpCensorIndicator!=0]
		tmpLower[censored>0] <- object$ucl[idx:iStart][tmpCensorIndicator>0]
		tmpUpper[censored<0] <- object$lcl[idx:iStart][tmpCensorIndicator<0]
		ysim <- tmvtnorm::rtmvnorm(1, mean=tmpMean, sigma=tmpVar, lower=tmpLower, upper=tmpUpper,algorithm="gibbs")
		y[idx] <- ysim[1]
  }
	m2 <- carx(y, object$x, rep(0,nObs), NULL, NULL, object$p,
	           prmtrAR=object$prmtrAR,
	           prmtrX = object$prmtrX,
	           sigma = object$sigma,
	           CI.compute=FALSE)
	#print(m2)
	rsdl <- numeric(nObs)
	eta <- y - object$x%*%m2$prmtrX
	for(idx in (p+1):nObs)
	{
		rsdl[idx] <- eta[idx] - eta[(idx-1):(idx-p)]%*%m2$prmtrAR
	}
	if(type=="raw")
		return(rsdl)
	if(type=="pearson")
	        return(rsdl/m2$sigma)
}

