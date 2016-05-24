#' Check the resolution of the data
#'
#' @description The function provides the resolution of the dendrometer data.
#'
#' @usage dendro.resolution(dm.data, unts = c("secs", "mins", "hours", "days"))
#'
#' @param dm.data a \code{data.frame} with a timestamp (\code{\%Y-\%m-\%d \%H:\%M:\%S} format) as row names, and dendrometer series in columns. Output as created using code from the \code{Import dendrometer data} vignette.
#' @param unts a \code{character} string of "secs", "mins", "hours", "days", specifiying the units in which the resolution should be calculated. Defaults to \code{"secs"}. Argument matching is performed.
#'
#' @return The function returns the resolution of the data in the desired unit.
#'
#' @author Marko Smiljanic
#'
#' @examples
#' data(dmCD, dmHS, dmED)
#' dendro.resolution(dmCD, unts = "hours")
#' dendro.resolution(dmHS, unts = "hours")
#' dendro.resolution(dmED, unts = "mins")
#'
#' @export
#'
dendro.resolution <- function(dm.data, unts = c("secs", "mins", "hours", "days"))
{
  nm <- deparse(substitute(dm.data))
  if(!is.dendro(dm.data)) {
		stop(paste("'", nm, "' is not in the required format", sep = ""))
  }

  unts <- match.arg(unts, c("secs", "mins", "hours", "days"))
	rowNames <- row.names(dm.data)
	rowdiff <- diff(as.POSIXct(rowNames, tz = "GMT"), units = unts)
	if(units(rowdiff) != unts) {
		units(rowdiff) <- unts
	}
	resolution <- unique(rowdiff)
	return(resolution)
}
