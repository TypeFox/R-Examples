#' Confirm that the metric functions are in suitable format
#'
#' Utility function. Creates a list of functions, either those defined in defineMetrics
#' or a named list of metric functions.
#'
#' @param x Optional named list of metric functions. Else, defines the metrics as those
#' defined in defineMetrics.
#' 
#' @details A few quick checks to confirm the metric functions are input in suitable
#' format.
#'
#' @return A list of functions.
#'
#' @export
#'
#' @references Miller, E. T., D. R. Farine, and C. H. Trisos. 2015. Phylogenetic community
#' structure metrics and null models: a review with new methods and software.
#' bioRxiv 025726.
#'
#' @examples
#' checkBetaMetrics(defineBetaMetrics())

checkBetaMetrics <- function(x)
{
	if (is.null(x))
	{
		metrics <- defineBetaMetrics()
	}
	else
	{
		if (!inherits(x, "list"))
		{
			stop("The metrics need to be input as a list of named functions")
		}
		if (is.null(names(x)))
		{
			stop("The metrics need to be input as a list of named functions")
		}
				
		metrics <- x
	}
	
	metrics
}
