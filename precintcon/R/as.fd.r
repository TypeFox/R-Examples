#' @include precintcon.fd.r
#' @include as.precintcon.fd.r
NULL

#' @name as.fd
#' @author Lucas Venezian Povoa \email{lucasvenez@@gmail.com} 
#' @aliases as.precintcon.fd as.fd precintcon.fd 
#' @title Frequency distribution of a precipitation serie 
#' @description Calculates the frequency distribution of a daily precipitation 
#' serie based on \code{interval}.
#' @usage as.precintcon.fd(object, interval = 1) 
#' @param object a daily precipitation serie.
#' @param interval the interval in millimeters for calculating the frequency 
#' distribution.
#' @return A data.frame (precintcon.fd) containing the following variables:
#' \itemize{
#' \item \code{initial.class} is the initial value of the class.
#' \item \code{final.class} is the final value of the class.
#' \item \code{midpoint} is the middle point of the class.
#' \item \code{n} is the absolute freqnecy, i.e., the number of days in each 
#' class.
#' \item \code{sum.n} is the cumulative frequency, obtained by adding the absolute 
#' frequencies of all the classes up to the one under consideration.
#' \item \code{P} is the pluviometric total of each class, obtained by multiplying 
#' \code{midpoint} by \code{n}.
#' \item \code{sum.P} is the cummulative class's pluviometric total, obtained by 
#' adding the pluviometric total of all the classes up to the one under consideration.
#' \item \code{p.sum.n} is the cumulative percentage of rainy days.
#' \item \code{p.sum.P} is the cumulative percentage of rainfall amounts.
#' }
#' @seealso
#' \code{\link{read.data}}
#' \code{\link{as.daily}}
#' \code{\link{ci}}
#' @keywords precipitation frequency distribution
#' @examples 
#' ##
#' # Loading the daily precipitation serie
#' data(daily)
#' 
#' ##
#' # Performing the frequency distribution
#' as.precintcon.fd(daily)
#' @export 
as.fd <- function(object, interval = 1) {
   return(as.precintcon.fd(object, interval))   
}