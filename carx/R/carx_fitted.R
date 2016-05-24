#' Fitted values of a \code{carx} object
#'
#' Compute the fitted values from a \code{carx} object.
#'  Note that the existence of censoring invalidates the usual Markov property for an AR model.  
#'  Instead, the conditional distribution of \eqn{Y^*_t} given the past \eqn{Y}s and current
#'  and past covariates is the same as the conditional distribution
#'  \eqn{D_t = D(Y^*_t|X_t, {(Y_{j}, X_j )}_{j=\tau}^{t-1} )},
#' where \eqn{1 \le \tau \le t-p} is the largest integer \eqn{t} such that 
#' none of \eqn{Y_t;t=\tau+p-1,...,\tau} is censored. In the case that \eqn{\tau = t-p}, 
#' the fitted value can be readily computed; otherwise, the fitted value is computed as the 
#' mean of the distribution 
#' \eqn{D_t} by the function \code{mtmvnorm} from the package \pkg{{tmvtnorm}}.

#' @param object a fitted \code{carx} object.
#' @param ... not used.
#' @return the fitted values.
#' @export
#' @examples
#' dat = carxSim(nObs=100,seed=0)
#' mdl <- carx(y~X1+X2-1,data=dat, p=2, CI.compute = FALSE)
#' #compute the fitted values
#' fv = fitted(mdl)
#
fitted.carx <- function(object,...)
{
	#message("Calling fitted.carx")
	nObs  <- object$nObs
	p  <- object$p

	ret <- rep(NA,nObs)
	trend <- object$x%*%object$prmtrX
	eta <- object$y - trend

	for(idx in (1:nObs)[-object$skipIndex])
	{
		#message(sprintf("calculating %i",idx))
		if(all(object$ci[(idx-1):(idx-p)]==0))
		{
			ret[idx] <- trend[idx] + eta[(idx-1):(idx-p)]%*%object$prmtrAR
		}
		else
		{
			#message(sprintf("Index %i is censored",idx))
			#find the beginning index of p consecutive observations
			iStart <- 1
			for(i in (idx-p):1)
			{
				if(all(object$ci[i:(i+p-1)]==0))
				{
					iStart <- i
					break
				}
			}

			nStart <- idx - iStart
			#message(sprintf("idx: %i, iStart: %i, nStart: %i",idx,iStart,nStart))
			tmpCensorIndicator <- object$ci[(idx-1):iStart] #looking back
			nCensored <- sum(tmpCensorIndicator!=0)
			covEta <- computeCovAR(object$prmtrAR, object$sigma, nStart+1)
			if( nCensored < nStart )
			{
				conditionalIndex <- which(tmpCensorIndicator==0) + 1
				tmpY <- object$y[idx:iStart][conditionalIndex]
				cdist <- conditionalDistMvnorm(tmpY, conditionalIndex, trend[idx:iStart], covEta)
				tmpMean <- cdist$'mean'
				tmpVar <- cdist$'var'
			}else
			{
				tmpMean <- trend[idx:iStart]
				tmpVar <- covEta
			}

			tmpLower <- rep(-Inf,length = nCensored+1) #( y[idx], censored obs)
			tmpUpper <- rep(Inf,length = nCensored+1)
			censored <- tmpCensorIndicator[tmpCensorIndicator!=0]
			tmpLower[-1][censored>0] <- object$ucl[(idx-1):iStart][tmpCensorIndicator>0]
			tmpUpper[-1][censored<0] <- object$lcl[(idx-1):iStart][tmpCensorIndicator<0]
			mn <- tmvtnorm::mtmvnorm(mean=tmpMean, sigma=tmpVar,lower=tmpLower,upper=tmpUpper,doComputeVariance=FALSE)
			ret[idx] <- mn$tmean[1]
		}
	}
	ret
}
