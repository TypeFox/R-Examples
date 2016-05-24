#' qm
#' 
#' Wraps an input in quotation marks.
#' 
#' @param x Character or Numeric.
#' 
#' @return A character string that begins and ends with quotation marks.
#' @author Jonathan A. Greenberg
#' 
#' @examples {
#' qm("Hi!")
#' qm(42) 
#' }
#' @export

qm <- function(x)
{
	paste('"',x,'"',sep="")
}