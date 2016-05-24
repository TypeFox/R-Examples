#' Simulate  data from a \code{carx} model
#'
#' Use the provided parameters in the supplied \code{carx} model and other settings to 
#' simulate data from the \code{carx} model; see  Wang and Chan (2015). 

#' @references Wang C, Chan KS (2015). "Quasi-likelihood estimation of a censored autoregressive model 
#' with exogenous variables." Submitted.
#'
#' @seealso \code{\link{carx}} for model specification. 
#'
#' @param nObs number of observations to be simulated.
#' @param prmtrAR the AR parameter.
#' @param prmtrX the regression parameters for X.
#' @param sigma the innovation standard deviation for the AR process.
#' @param lcl the lower censoring limit.
#' @param ucl the upper censoring limit.
#' @param x optional matrix for X. Default = \code{NULL}, in which case X will be simulated from 
#'  the standard normal distribution with dimensions determined by \code{nObs} and \code{prmtrX}.
#' @param seed optional to set the seed of random number generator used by \code{R}, default=\code{NULL}.
#' @param inno.dist innovation distribution, can be "normal" or "t", default="normal". If it is "t", 
#'  its degree of freedom should be supplied in \code{t.df}.
#' @param t.df the degree of freedom of the t distribution, used only if \code{inno.dist}="t". Default=5.
#' @return a data frame of simulated \code{y}, \code{x}, \code{ci}, \code{lcl} and \code{ucl}.
#' @export
#' @examples
#' dat = carxSim()

carxSim <- function(nObs=200, prmtrAR=c(-0.28,0.25), prmtrX=c(0.2,0.4), sigma=0.60, lcl=-1, ucl=1, x = NULL, seed=NULL,inno.dist=c("normal","t"),t.df=5)
{
  if(is.null(lcl)) lcl = -Inf
  if(is.null(ucl)) ucl = Inf

  stopifnot(t.df>2)
	p <- length(prmtrAR)
	nX <- length(prmtrX)
  inno.dist = match.arg(inno.dist)

	if(!is.null(seed))
		set.seed(seed)

  if(inno.dist == "normal")
    eps <- stats::rnorm(nObs,0, sigma)
  else
  {
    if(inno.dist == "t")
    {
      eps <- stats::rt(nObs,t.df)
      eps <- eps*sigma/(sqrt(t.df/(t.df-2)))
    }
  }


	if(is.null(x))
  {
		x <- matrix(stats::rnorm(nObs*nX), nrow= nObs, ncol = nX)
  }

  if(is.null(colnames(x)))
    colnames(x) <- paste0("X",seq(1,nX))

	trend <- x%*%prmtrX

	eta <- numeric(nObs)
	y <- numeric(nObs)


  if(inno.dist == "normal")
  {
    covAr <- computeCovAR(prmtrAR,sigma,lag=p)
    eta[1:p] <- as.vector(mvtnorm::rmvnorm(1,sigma=covAr))
  }
  else
  {
    if(inno.dist == "t")
    {
      nPreSample <- 1000
      tmpEps <- stats::rt(nPreSample,t.df)
      tmpEps <- tmpEps*sigma/(sqrt(t.df/(t.df-2)))
      tmpEta <- rep(0,nPreSample)
      for(i in (p+1):nPreSample)
        tmpEta[i] <- prmtrAR%*%tmpEta[(i-1):(i-p)] + tmpEps[i]
      eta[1:p] <- tmpEta[(nPreSample-p+1):nPreSample]
    }
  }


  #assign first p values
	y[1:p] <- trend[1:p] + eta[1:p]

  #iterate
	for(i in (p+1):nObs){
		eta[i] <- eta[(i-1):(i-p)] %*% prmtrAR + eps[i]
		y[i] <- trend[i] + eta[i]
	}
	ci <- rep(0,nObs)
	ci[y<lcl] <- -1
	ci[y>ucl] <- 1

  if(options()$verbose) message(paste0("simulated series: censor rate: ", sum(abs(ci))/nObs))
	ret <- list(y = y,
		    ci=ci,
		    lcl=lcl,
		    ucl=ucl,
		    x = x
		    )
  ret <- try(data.frame(ret,row.names=NULL),silent=TRUE)
  #assign("last.warning", NULL, envir = baseenv())
  colnames(ret) <- c("y","ci","lcl","ucl",colnames(x))
	ret
}


#' simulate a sample \code{\link{cenTS}} data for \code{carx}
#'
#' Use provided parameters and other settings to simulate a series of data as a \code{cenTS} object.
#'
#' @inheritParams carxSim
#' @param value.name the name of the response series
#' @param end.date the date of the last observation, default = \code{Sys.date()}.
#' @return a \code{cenTS} object with regressors.
#' @seealso \code{\link{carxSim}}.
#' @export
#' @examples
#' cts = carxSimCenTS()
carxSimCenTS <- function(nObs=200, prmtrAR=c(-0.28,0.25), prmtrX=c(0.2,0.4), sigma=0.60, lcl=-1, ucl=1, x = NULL, seed=NULL, value.name = 'y', end.date=Sys.Date(),inno.dist=c("normal","t"),t.df=5)
{
    ret <- carxSim(nObs,prmtrAR, prmtrX, sigma, lcl, ucl, x, seed)
  #ret is a data.frame
  names(ret) <- c("value",names(ret)[-1])
  ret <- as.list(ret)
  ret$order.by <- end.date+seq(-nObs,-1,by=1)
  #listx$value <- ret$y
  ret$value.name <- value.name

  val <- do.call(cenTS,ret)
  val
}
