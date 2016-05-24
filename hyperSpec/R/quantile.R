.quantile <- function (x, probs = seq(0, 1, 0.5), na.rm = TRUE, names = "num", ...){
	
	x <- apply (x, 2, quantile, probs = probs, na.rm = na.rm, names = FALSE, ...,
							long = list (probs = probs, na.rm = na.rm, names = names, ...))
	
	if (names == "pretty") 
		rownames (x@data) <- paste (format (100 * probs, format = "fg", width = 1,
																				justify = "right",
																				digits =  getOption ("digits")),
																"%")
	else if (names == "num")
		rownames (x@data) <- probs
	
	x
}

##' @rdname mean_sd
##' @return For hyperSpec object, \code{quantile} returns a hyperSpec object containing the
##' respective quantile spectra.
##' @param probs the quantiles, see \code{\link[stats]{quantile}}
##' @param names \code{"pretty"} results in percentages (like \code{\link[stats]{quantile}}'s
##' \code{names = TRUE}), \code{"num"} results in the row names being \code{as.character (probs)}
##' (good for ggplot2 getting the order of the quantiles right). Otherwise, no names are assigned.
##' @seealso  \code{\link[stats]{quantile}}
##' @export
##' @examples
##' 
##' plot (quantile (chondro))
setMethod ("quantile", signature = signature (x = "hyperSpec"), .quantile)
