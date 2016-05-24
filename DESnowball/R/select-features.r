#' Compute ranking statistics, RD, and p value for gene selection
#'
#' Gene selection based on the statistical significances according to the Snowball approach (see references for more details).
#' @param x an output from the main function \code{\link{snowball}}
#' @param cutoff.p cutoff for top gene list. This is applied on the multiple testing adjusted p values
#' @param p.adjust.method specifies a multiple testing adjustment method, see \code{\link{p.adjust}} for more details
#' @return a list with two elements - \code{fullList} and \code{selectedList}.\code{fullLIst} 
#' is a data.frame that contains \code{rd}, \code{pval} and \code{positive}, corresponding respectively 
#' to the RD, p value and an indicator variable of weather the RD value is above or below the 
#' median value. \code{selectedList} is a data.frame that contains the same variables as those 
#' in \code{fullList} with only the top genes that satisfy the significance cutoff specified 
#' by \code{cutoff.p}.
#' @references
#' Xu, Y., Guo, X., Sun, J. and Zhao. Z. Snowball: resampling combined with distance-based regression to discover transcriptional consequences of driver mutation, manuscript.
#' @export
select.features <- function(x,
			    cutoff.p=0.05,
			    p.adjust.method='BH')
{
    full.list <- fs.rd(x=x,
		       method="mcd",
		       df=1)

    selected.list <- fs.selection(full.list,
				  cutoff.p = cutoff.p,
				  p.adjust.method=p.adjust.method)
    ret <- list(fullList=full.list,
		selectedList=selected.list)
    class(ret) <- 'fsfeature'
    ret
}
