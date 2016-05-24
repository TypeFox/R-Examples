#' Simulate a significance threshold for (composite) interval mapping
#' 
#' Generate a significance threshold by simulating from the null hypothesis. Phenotypic values for the observed genetic data are simulated to represent data with no QTL, and a genomewide p-value threshold is calculated. 
#' @export
#' @param mpcross Object of class \code{mpcross}
#' @param nsim Number of null datasets to simulate; default is 100
#' @param alpha Significance threshold for QTL p-values
#' @param pindex Index of phenotype (if more than one in dataset)
#' @param step Step size at which to compute the QTL profile. Default is 0, at midpoints of marker intervals
#' @param ncov Number of marker covariates to search for - default is 0 for interval mapping
#' @param ... Additional arguments
#' @return List with components:
#' \item{alpha}{Input significance threshold}
#' \item{nsim}{Input number of simulated datasets}
#' \item{minp}{Genomewide minimum p-value for each simulated dataset}
#' \item{thr}{Empirical p-value threshold}
#' @seealso \code{\link[mpMap]{mpIM}}

sim.sigthr <-
function(mpcross, nsim=100, alpha=.05, pindex=1, step=0, ncov=0, ...)
{
  output <- list()
 
  if (is.null(mpcross$map)) 
	stop("Must have marker map to generate null simulations")

  vare <- 1
  if (!is.null(mpcross$pheno))
	vare <- var(mpcross$pheno[,pindex])

	#if the ibd probabilities haven't been computed yet, compute them
	if (!(inherits(mpcross, "mpprob") && attr(mpcross$prob, "step")==step))
		mpp <- mpprob(mpcross, program="qtl", step=step)
	else mpp <- mpcross

  minp <- vector(length=nsim)
  for (i in 1:nsim)
  {
#    dat <- sim.mpcross(mpcross$map, mpcross$ped, vare=vare, seed=i)
    ph <- rnorm(nrow(mpcross$finals))
    mpp$pheno[,1] <- ph
    res <- mpIM(object=mpp, ncov=ncov, step=step, responsename="pheno", ...)
    minp[i] <- min(unlist(res$QTLresults$pvalue))
  }
 
  output$call <- match.call()
  output$minp <- minp
  output$thr <- sort(minp)[floor(alpha*nsim)]

  return(output)
}

