#' @title Compute log-likelihood of a meteDist object
#'
#' @description
#' \code{logLik.meteDist} computes log-likelihood of a meteDist object
#'
#' @details
#' Degrees of freedom are assumed to be equal to the number of Lagrange 
#' multpliers needed to specify the METE prediction. See Examples for usage.
#' 
#' 
#' @param object a \code{meteDist} object
#' @param ... arguments to be passed
#' @export
#' 
#' @examples
#' data(arth)
#' ## object holding ecosystem structure function
#' esf1 <- meteESF(spp=arth$spp,
#'                 abund=arth$count,
#'                 power=arth$mass^(.75),
#'                 minE=min(arth$mass^(.75)))
#' ## calculate individual power distribution and its likelihood
#' ipd1 <- ipd(esf1)
#' logLik(ipd1)
#' 
#' @return object of class \code{logLik}
#'
#' @author Andy Rominger <ajrominger@@gmail.com>, Cory Merow
#' @seealso sad, ssad, ipd, sipd
#' @references Harte, J. 2011. Maximum entropy and ecology: 
#' a theory of abundance, distribution, and energetics. Oxford University Press.

logLik.meteDist <- function(object,...) {
  lik <- sum(object$d(object$data,log=TRUE))
  attr(lik, 'df') <- length(object$La)
  class(lik) <- 'logLik'
  return(lik)
}

#==============================================================================
#' @title Compute log-likelihood z-score
#'
#' @description
#' \code{logLikZ.meteDist} computes a log-likelihood z-score by simulation from a 
#' fitted METE distribution 
#'
#' @details
#' \code{logLikZ.meteDist} simulates from a fitted METE distribution (e.g. a species 
#' abundance distribution or individual power distribution) and calculates the 
#' likelihood of these simulated data sets. The distribution of these values is compared 
#' against the likelihood of the data to obtain a z-score. 
#' 
#' @param x a \code{meteDist} object
#' @param nrep number of simulations from the fitted METE distribution 
#' @param return.sim logical; return the simulated liklihood values
#' @param ... arguments to be passed to methods
#' @export
#' 
#' @examples
#' data(arth)
#' ## object holding ecosystem structure function
#' esf1 <- meteESF(spp=arth$spp,
#'                 abund=arth$count,
#'                 power=arth$mass^(.75),
#'                 minE=min(arth$mass^(.75)))
#' ## calculate individual power distribution
#' ipd1 <- ipd(esf1)
#' ## calculate z-score, keeping all simulated log likelihoods for plotting
#' llz <- logLikZ(ipd1, nrep=100, return.sim=TRUE)
#' 
#' plot(density(llz$sim),xlim=range(c(llz$sim,llz$obs)),
#'      xlab='log(likelihood)',col='red')
#' abline(v=llz$obs,lty=2)
#' legend('top',legend=c('data','simulated'),col=c('black','red'),
#'       lty=c(1,1),bty='n') 
#' 
#' @return list with elements
#' \describe{
#'    \item{z}{The z-score}
#'    \item{obs}{The observed log-likelihood of the actual data set}
#'    \item{sim}{\code{nrep} Simulated values if return.sim=TRUE, NULL otherwise}
#' }
#'
#' @author Andy Rominger <ajrominger@@gmail.com>, Cory Merow
#' @seealso mseZ.meteDist
#' @references Harte, J. 2011. Maximum entropy and ecology: a theory of abundance, distribution, and energetics. Oxford University Press.

logLikZ <- function(x, ...) {
  UseMethod('logLikZ')
}

#' @rdname logLikZ
#' @export 

#' @importFrom stats logLik sd
logLikZ.meteDist <- function(x, nrep, return.sim=FALSE,...) {
  lik.obs <- logLik(x)
  lik.sim <- replicate(nrep, {
    new.dat <- x$r(length(x$data))
    sum(x$d(new.dat, log=TRUE))
  })
  
  z <- (lik.obs-mean(lik.sim))/sd(lik.sim)
  
  if(!return.sim) lik.sim <- NULL
  
  return(list(z=z, 
              obs=lik.obs,
              sim=lik.sim))
  
}
