


#' Custom rounding function to round to the nearest specified interval
#' @param x numeric value(s)
#' @param roundTo rounding interval
#' @return rounded numeric value(s)
#' @author Jason Waddell
#' @export
customRound <- function(x, roundTo){
	round(x/roundTo)*roundTo
}
