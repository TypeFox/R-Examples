#' @title 
#' Sample Size Search for EpiBayes Models
#' 
#' @description
#' This function takes three vectors for the number of subzones, clusters per subzone, and
#'     subjects per cluster per subzone and uses all combinations of the given values 
#'     as inputs to \code{\link{EpiBayes_ns}} as a way to search for sample sizes which give
#'     optimal values of the Bayesian model output (e.g., \code{p4.tilde}). It assumes that
#'     at every sampling level, all elements have the same size (e.g., if the user supplies
#'     H = 2, k = 10, n = 100, then we assume that there are two subzones, both of which
#'     contain 10 clusters/farms/ponds/herds, and all clusters in all subzones contain 
#'     100 subjects/individuals/mollusks/cows/chickens). This is done for computational
#'     efficiency in the search.
#'
#' @param H.vect Values of possible numbers of subzones. Integer vector.
#' @param k.vect Values of possible numbers of clusters within subzones. Integer vector.
#' @param n.vect Values of possible numbers of subjects within clusters within subzones. 
#'     Integer vector.
#' @param season.vect The single season in which one assumes sampling is taking place. Coded as 
#'     (1) Summer, (2) Fall, (3) Winter, (4) Spring. Integer scalar.
#' @param ... Additional arguments that will be passed to \code{\link{EpiBayes_ns}}. 
#'     Otherwise, the default values will be used.
#'
#' @return
#' The returned values are given in a matrix. They are as follows.
#'
#' \tabular{lll}{
#'     Output \tab Attributes \tab Description \cr
#'     \code{RawPost} \tab List: Length - (number of periods), Elements - Real arrays (\code{reps} x \code{H} x \code{MCMCreps}) \tab Posterior distributions for the cluster-level prevalences for each subzone from all time periods \cr
#'     \code{BetaBusterEst} \tab List: Length - (number of periods), Elements - Real vectors (2 x 1) \tab Estimated posterior distributions for the cluster-level prevalences for each subzone from all time periods using moment-matching to the closest beta distribution by the function \code{\link[epiR]{epi.betabuster}} \cr
#'     \code{ForOthers} \tab  \tab Various other data not intended to be used by the user, but used to pass information on to the \code{plot}, \code{summary}, and \code{print} methods \cr
#' }
#'
#' @examples
#' testrun_samplesize = EpiBayesSampleSize(
#'		H = c(2, 4), 
#'		k = c(10, 20), 
#'		n = c(100, 500),
#'		season = 3,
#'		burnin = 1,
#'		reps = 1,
#'		MCMCreps = 10,
#'		tau.T = 0,
#'		poi = "tau",
#'		mumodes = matrix(c(
#'			0.50, 0.70, 
#'			0.50, 0.70, 
#'			0.02, 0.50, 
#'			0.02, 0.50
#'			), 4, 2, byrow = TRUE
#'		),
#'		pi.thresh = 0.05, 
#'	    tau.thresh = 0.02,
#'      gam.thresh = 0.10,
#'		poi.lb = 0.1,
#'		poi.ub = 0.4,
#'		p1 = 0.95,
#'		psi = 4,
#'		tauparm = c(1, 1),
#'		omegaparm = c(1000, 1),
#'		gamparm = c(1000, 1),
#'		etaparm = c(100, 6),
#'		thetaparm = c(100, 6) 
#'		)
#'
#' testrun_samplesize
#' print(testrun_samplesize, out.ptilde = "p4.tilde")
#' 
#' @export
#' @name EpiBayesSampleSize
#' @import compiler
#' @import epiR 
EpiBayesSampleSize = function(
						H.vect,
						k.vect,
						n.vect,
						season.vect,
			      		...
						){
	
	## Create a grid of H/k/n values to compute over
	searchgrid = expand.grid(H.vect, k.vect, n.vect)
	H.list = list()
	k.list = list()
	n.list = list()
	season.list = list()

	for(i in 1:nrow(searchgrid)){
		H.list[[i]] = searchgrid[i, 1]
		k.list[[i]] = rep(searchgrid[i, 2], searchgrid[i, 1])
		n.list[[i]] = rep(rep(searchgrid[i, 3], searchgrid[i, 2]), searchgrid[i, 1])
		season.list[[i]] = rep(rep(season.vect, searchgrid[i, 2]), searchgrid[i, 1])
	}

	
	## Store the extra arguments to EpiBayes_ns or EpiBayes_s into a variable to pass to mapply()
	EBargs = list(...)
	
	# Store the lower and upper bounds for the print method
		if (!is.null(EBargs$poi.lb)){
			poi.lb = EBargs$poi.lb
		} else{
			poi.lb = 0
			}
			
		if (!is.null(EBargs$poi.ub)){
			poi.ub = EBargs$poi.ub
		} else{
			poi.ub = 1
			}
	
	## Use mapply to search over the grid
	samplesearch = mapply(
						EpiBayes_ns,
						H = H.list,
						k = k.list,
						n = n.list,
						seasons = season.list,
						MoreArgs = EBargs
						)							
						
	EpiBayesSampleSize.out = list(
								"samplesearch" = samplesearch,
								"ForOthers" = list(
												"H.vect" = H.vect,
												"k.vect" = k.vect,
												"n.vect" = n.vect,
												"season.vect" = season.vect,
												"searchgrid" = searchgrid,
												"poi.lb" = poi.lb,
												"poi.ub" = poi.ub
												))
								
	## Store the object in its own class so it has a print method at least
	class(EpiBayesSampleSize.out) = "ebsamplesize"
	
	return(invisible(EpiBayesSampleSize.out))
}