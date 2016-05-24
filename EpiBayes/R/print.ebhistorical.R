#' @title 
#' Print Method for EpiBayes Historical Object
#' 
#' @description
#' This function gives basic summary statistics for posterior distributions of cluster-level 
#'     prevalences across all time periods considered. It does so by examining the object 
#'     output by the \code{\link{EpiBayesHistorical}} function of class
#'     \code{ebhistorical}.
#'
#' @param x An object of class \code{ebhistorical} (e.g., the output of 
#'     function \code{\link{EpiBayesHistorical}}).
#' @param ... Additional arguments to be passed on to summary.
#'
#' @return
#' The posterior means are returned in a matrix with rows indexed by time period and 
#'     columns indexed by subzone. The output of this method can be customized using
#'     the \code{summary} method for this \code{ebhistorical} class of object.
#'
#' @seealso 
#' This is a method for objects of class \code{ebhistorical} returned by the function
#'     \code{\link{EpiBayesHistorical}}. The \code{summary} method for
#'     the \code{ebhistorical} object class allows for customization of this output.
#'
#' @export
#' @name print.ebhistorical
print.ebhistorical = function(x, ...){
		
		# Extract relevant information
		periodi.tau = x$RawPost
		periodi.betabust = x$BetaBusterEst
		unique.subzones = x$ForOthers$unique.subzones
		burnin = x$ForOthers$burnin
		
		# Get dimensions for output
		n.periods = length(periodi.betabust)
		n.subzones = nrow(periodi.betabust[[1]])
		
		# Compute summary statistic for each subzone	
		sumstat.mat = matrix(unlist(lapply(periodi.tau,
					   function(x){
					       apply(x, 
					       		 1, 
					       		 function(y){
					       		 			mean(y[-c(1:burnin)])
					       		 }
					       	)
					    }
					)
				),
				n.subzones, n.periods)			
		
		sumstat.mat = list(
						"sumstat.mat" = sumstat.mat, 
						"sumstat.call" = match.call(expand.dots = TRUE), 
						"ForOthers" = list(
										unique.subzones = unique.subzones
										)
						)
		
		return(sumstat.mat)
							
	}