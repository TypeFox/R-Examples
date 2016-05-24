#' Check input data
#'
#' @description The function checks whether the input data is in the required format, as described in the \code{Import} \code{dendrometer} \code{data} vignette.
#'
#' @usage is.dendro(dm.data)
#'
#' @param dm.data a \code{data.frame} with a timestamp (\code{\%Y-\%m-\%d \%H:\%M:\%S} format) as row names, and dendrometer series in columns. Output as created using code from the \code{Import dendrometer data} vignette.
#'
#' @return The function returns \code{TRUE} if the input data is valid and \code{FALSE} otherwise. In the latter case, specific error messages are given as well.
#'
#' @author Ernst van der Maaten, Marieke van der Maaten-Theunissen and Marko Smiljanic.
#'
#' @examples
#' data(dmCD, dmHS, dmED)
#' is.dendro(dmCD)
#' is.dendro(dmHS)
#' is.dendro(dmED)
#'
#' @importFrom methods is
#' @export
#'
is.dendro <- function(dm.data)
{
  nm <- deparse(substitute(dm.data))

	if(!is.data.frame(dm.data)) {
		warning(paste("'", nm, "' is not a data frame", sep = ""))
	  return(FALSE)
	}
  tst.timestamp.data <- try(is(as.POSIXct(rownames(dm.data)),"POSIXct"), silent = TRUE)
  if(tst.timestamp.data != TRUE){
    warning(paste("rownames of '", nm, "' is not a timestamp or contains errors", sep = ""))
    return(FALSE)
  }
  if(length(rownames(dm.data)) != length(unique(rownames(dm.data)))) {
    warning(paste("the date-time stamp of '", nm, "' contains non-unique values", sep = ""))
    return(FALSE)
  }
  row.nm <- row.names(dm.data)
  row.diff <- diff(as.POSIXct(row.nm, tz = "GMT"))
  if(units(row.diff) != "hours") {
    units(row.diff) <- "hours"
  }
  resolution <- unique(row.diff)
  if(length(resolution) != 1) {
    warning(paste("the temporal resolution of '", nm, "' is not constant", sep = ""))
    return(FALSE)
  }
  tst.data.numeric <- sapply(dm.data, is.numeric)
  if(FALSE %in% tst.data.numeric) {
    warning(paste("the columns of '", nm, "' should contain numeric dendrometer data only", sep = ""))
    return(FALSE)
  }
	return(TRUE)
}
