#' @title 
#' Print Method for EpiBayes Historical Object
#' 
#' @description
#' This function prints the output of objects of class \code{ebhistoricalsummary} from the 
#'     method \code{\link{summary.ebhistorical}} in a nice form.
#'
#' @param x An object of class \code{ebhistoricalsummary} (e.g., the output of 
#'     function \code{\link{summary.ebhistorical}}).
#' @param ... Additional arguments to be passed on to print.
#'
#' @return
#' The summary statistics are returned in a matrix with rows indexed by time period and 
#'     columns indexed by subzone. The body of the matrix is filled with the summary 
#'     statistics requested for by the argument \code{sumstat}.
#'
#' @seealso 
#' This is a method for objects of class \code{ebhistoricalsummary} returned by the method
#'     \code{\link{summary.ebhistorical}}, which creates its own class of object much like 
#'     the summary method for \code{lm} objects. There is an additional plot method that 
#'     will plot the summary output called \code{\link{plot.ebhistoricalsummary}}.
#'
#' @export
#' @name print.ebhistoricalsummary
print.ebhistoricalsummary = function(x, ...){
	
	# Get some labels from the historicalsummary object
	sumstat.mat = x$sumstat.mat
	sumstat.call = x$sumstat.call
	unique.subzones = x$ForOthers$unique.subzones
	time.labels = x$ForOthers$time.labels
	n.periods = ncol(sumstat.mat)
	
		# Keep the subzone names			
		rownames(sumstat.mat) = unique.subzones
		
		# Keep the time period names
		if (is.null(time.labels)){
			time.labels = 1:n.periods
		}	
		
		colnames(sumstat.mat) = time.labels
		
	# Handle differently if the user has requested posterior quantiles
	# if (sumstat.call$sumstat == "quantile"){
		# cat(paste0("Matrix of posterior ", sumstat.call$prob, "th", sumstat.call$sumstat, "s of cluster-level prevalence across ", nrow(sumstat.mat), " time periods \n \n"))
	# } else{
		cat(paste0("Matrix of posterior ", sumstat.call$sumstat, "s of cluster-level prevalence across ", n.periods, " time periods \n \n"))
		# }
		
	print(sumstat.mat)
	
	return(sumstat.mat)

}