#' @name tii
#' @aliases tii
#' @author Lucas Venezian Povoa
#' 
#' @title Temporaly Irregularity Index
#' @description Calculates the temporaly irregularity index according to the function sum(log(Pi+1/Pi))/(n-1),
#' where Pi is the precipitation amount of year i, and n is the number of years.
#' @details Daily or monthly precipitation series are transformed to annual series.
#' @usage tii(object)
#' @param object is a daily or monthly precipitation serie
#' @return the temporaly irregularity index according to the function sum(log(Pi+1/Pi))/(n-1)
#' @examples
#' ##
#' # Loading the monthly precipitation serie
#' data(monthly)
#' 
#' ##
#' # Calculationg the Temporaly Irregularity Index
#' tii(monthly)
#' @export
tii <- function(object) {
  object <- as.annual(object)
  return(sum(log(object$precipitation[2:nrow(object)]/object$precipitation[(1:(nrow(object)-1))]))/(nrow(object)-1))
}
