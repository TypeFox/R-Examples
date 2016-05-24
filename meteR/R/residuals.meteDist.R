#' @title Compute residuals between METE predictions and date of a meteDist object
#'
#' @description
#' \code{residuals.meteDist} computes residuals between METE predictions and 
#' data of a meteDist object
#'
#' @details
#' See Examples. Typically not called directly by the user and rather used for 
#' calculating the mean square error with \code{mse.meteDist}. If \code{type='rank'}
#' returned value will be of length equal to number of observations (e.g. number of
#' species in case of SAD) but if \code{type='cumulative'} returned value will be of
#' length equal to number of unique ovservations (e.g. number of unique abundances in
#' case of SAR).
#' 
#' @param object a \code{meteDist} object
#' @param type 'rank' or 'cumulative'
#' @param relative logical; if true use relative MSE
#' @param log logical; if TRUE calculate MSE on logged distirbution. If FALSE use arithmetic scale.
#' @param ... arguments to be passed to methods
# @keywords manip
#' @export
#' 
#' @examples
#' data(arth)
#' esf1 <- meteESF(spp=arth$spp,
#'                 abund=arth$count,
#'                 power=arth$mass^(.75),
#'                 minE=min(arth$mass^(.75)))
#' sad1 <- sad(esf1)
#' residuals(sad1)
#' 
#' @return a numeic vector giving residuals for each data point
#'
#' @author Andy Rominger <ajrominger@@gmail.com>, Cory Merow
#' @seealso \code{mse.meteDist}
#' @references Harte, J. 2011. Maximum entropy and ecology: a theory of abundance, distribution, and energetics. Oxford University Press.
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the documentation.

residuals.meteDist <- function(object, type=c("rank","cumulative"),
                               relative=TRUE, log=FALSE, ...) {
	type <- match.arg(type, choices=c("rank","cumulative"))
	
	if(type=="rank") {
		obs <- object$data
		pred <- meteDist2Rank(object)
	} else if(type=="cumulative") {
		obs <- .ecdf(object$data)
		pred <- object$p(obs[,1], log.p=log)
		if(log) obs <- log(obs[,2])
		else obs <- obs[,2]
	}
	
	out <- obs - pred
	if(relative) out <- out/abs(pred)
	
	return(out)
} 

#==========================================================================
#' @title Computes mean squared error for rank or cdf
#'
#' @description
#' \code{mse.meteDist} computes mean squared error for rank or cdf between METE prediction and data
#'
#' @details
#' See Examples.
#' 
#' @param x a \code{meteDist} object
#' @param type 'rank' or 'cumulative'
#' @param relative logical; if true use relative MSE
#' @param log logical; if TRUE calculate MSE on logged distirbution. If FALSE use arithmetic scale.
#' @param ... arguments to be passed to methods
# @keywords manip
#' @export
#' 
#' @examples
#' data(arth)
#' esf1 <- meteESF(spp=arth$spp,
#'                 abund=arth$count,
#'                 power=arth$mass^(.75),
#'                 minE=min(arth$mass^(.75)))
#' sad1 <- sad(esf1)
#' mse(sad1, type='rank', relative=FALSE)
#' ebar1 <- ebar(esf1)
#' mse(ebar1)
#' 
#' @return numeric; the value of the mean squared error.
#'
#' @author Andy Rominger <ajrominger@@gmail.com>, Cory Merow
#' @seealso mseZ.meteDist
#' @references Harte, J. 2011. Maximum entropy and ecology: a theory of abundance, distribution, and energetics. Oxford University Press.
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the documentation.

mse <- function(x, ...) {
  UseMethod('mse')
}

#' @rdname mse
#' @export 
#' @importFrom stats residuals 
mse.meteDist <- function(x, type=c("rank","cumulative"),
                         relative=TRUE, log=FALSE, ...) {
  type <- match.arg(type, choices=c("rank","cumulative"))
  
  resid <- residuals(x, type, relative, log)
  
  return(mean(resid^2))
}

#' @export 
#' @importFrom stats residuals
mse.meteRelat <- function(x,...) {
  resid <- residuals(x)
  
  return(mean(resid^2))
}

#=======================================================================
#' @title Compute z-score of mean squared error
#'
#' @description
#' \code{mseZ.meteDist} Compute z-score of mean squared error
#'
#' @details
#' \code{mseZ.meteDist} simulates from a fitted METE distribution (e.g. a species abundance distribution or individual power distribution) and calculates the MSE between the simulated data sets and the METE prediction. The distribution of these values is compared against the MSE of the data to obtain a z-score. 
#' 
#' @param x a \code{meteDist} object
#' @param nrep number of simulations from the fitted METE distribution 
#' @param return.sim logical; return the simulated liklihood values
#' @param type either "rank" or "cumulative"
#' @param ... arguments to be passed to methods
# @keywords manip
#' @export
#' 
#' @examples
#' esf1=meteESF(spp=arth$spp,
#'               abund=arth$count,
#'               power=arth$mass^(4/3),
#'               minE=min(arth$mass^(4/3)))
#' sad1=sad(esf1)
#' mseZ(sad1, nrep=100)
#' @return list with elements
#' \describe{
#'    \item{z}{The z-score}
#'    \item{obs}{Mean squared error}
#'    \item{sim}{\code{nrep} Simulated values}
#' }
#'
#' @author Andy Rominger <ajrominger@@gmail.com>, Cory Merow
#' @seealso logLik.meteDist
#' @references Harte, J. 2011. Maximum entropy and ecology: a theory of abundance, distribution, and energetics. Oxford University Press.
# @aliases - a list of additional topic names that will be mapped to
# this documentation when the user looks them up from the command
# line.
# @family - a family name. All functions that have the same family tag will be linked in the documentation.

mseZ <- function(x, ...) {
  UseMethod('mseZ')
}

#' @rdname mseZ
#' @export 
#' @importFrom stats sd
mseZ.meteDist <- function(x, nrep, return.sim=FALSE,
                          type=c("rank","cumulative"), ...) {
  #seems like this should default to relative resid, since that's what is produced by residuals() for the data. need to add an arg for relative=TRUE?s
  if(type=='rank') {
    thr <- function(dat) {
      mean((sort(dat, TRUE) - meteDist2Rank(x))^2)
    }
  } else {
    thr <- function(dat) {
      obs <- .ecdf(dat)
      mean((obs - x$p(obs[,1]))^2)
    }
  }
  
  mse.obs <- mse.meteDist(x, type)
  mse.sim <- replicate(nrep, {
    new.dat <- x$r(length(x$data))
    
    thr(new.dat)
  })
  
  if(return.sim) {
    return(list(z=(mse.obs-mean(mse.sim))/sd(mse.sim), 
                obs=mse.obs,
                sim=mse.sim))
  } else {
    return(list(z=(mse.obs-mean(mse.sim))/sd(mse.sim)))
  }
}

