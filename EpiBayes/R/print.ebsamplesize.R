#' @title 
#' Print Method for EpiBayes Sample Size Object
#' 
#' @description
#' This function prints the output of objects of class \code{ebsamplesize} from the 
#'     method \code{\link{EpiBayesSampleSize}} in a nice form.
#'
#' @param x An object of class \code{ebsamplesize} (e.g., the output of 
#'     function \code{\link{EpiBayesSampleSize}}).
#' @param out.ptilde Indicators for the desired output of the printing. Can be one or any 
#"     combination of: \code{p2.tilde}, \code{p4.tilde}, and/or \code{p6.tilde}. If 1 
#'     output is desired, supply a character vector of length 1, if 2 are desired, supply 
#'     a character vector of length 2, and 3 for 3. Character vector.
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
#' @name print.ebsamplesize
print.ebsamplesize = function(x, out.ptilde = c("p2.tilde", "p4.tilde", "p6.tilde"), ...){
	
	# Get some labels from the ebsamplesize object
	n.searches = ncol(x$samplesearch)
	H.vect = x$ForOthers$H.vect	
	k.vect = x$ForOthers$k.vect	
	n.vect = x$ForOthers$n.vect
	searchgrid = x$ForOthers$searchgrid	
	n.H = length(H.vect)
	n.k = length(k.vect)
	n.n = length(n.vect)
	
	# Get bounds for reporting
	poi.thresh = x$ForOthers$poi.thresh
	poi.lb = x$ForOthers$poi.lb
	poi.ub = x$ForOthers$poi.ub
	
	# Set up empty array to hold results
	p2table = array(0, dim = c(n.k, n.n, n.H), dimnames = list("k" = k.vect, "n" = n.vect, "H" = H.vect))
	p4table = array(0, dim = c(n.k, n.n, n.H), dimnames = list("k" = k.vect, "n" = n.vect, "H" = H.vect))
	p6table = array(0, dim = c(n.k, n.n, n.H), dimnames = list("k" = k.vect, "n" = n.vect, "H" = H.vect))
	
	# Fill in the tables
	for (i in 1:nrow(searchgrid)){
		curr.kind = which(dimnames(p2table)$k == searchgrid[i, 2])
		curr.nind = which(dimnames(p2table)$n == searchgrid[i, 3])
		curr.Hind = which(dimnames(p2table)$H == searchgrid[i, 1])
		
				p2table[curr.kind, curr.nind, curr.Hind] = x$samplesearch[, i]$p2.tilde
				p4table[curr.kind, curr.nind, curr.Hind] = x$samplesearch[, i]$p4.tilde
				p6table[curr.kind, curr.nind, curr.Hind] = x$samplesearch[, i]$p6.tilde

	}
	
		cat(paste0("Estimated probabilities of not detecting cluster-level prevalence above", poi.thresh, " for given H, k, and n. \n \n"))
		print(p2table)
	
		cat(paste0("Estimated probabilities of detecting cluster-level prevalence above", poi.thresh, " for given H, k, and n. \n \n"))
		print(p4table)
	
		cat(paste0("Estimated probabilities of detecting cluster-level prevalence in the interval (", poi.lb, ", ", poi.ub, ") for given H, k, and n. \n \n"))
		print(p6table)

	ebsamplesize.printout = list(
					"p2table" = p2table,
					"p4table" = p4table,
					"p6table" = p6table
					)

	return(invisible(ebsamplesize.printout))

}