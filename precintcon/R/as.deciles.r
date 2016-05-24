#' @include precintcon.deciles.analysis.r
#' @include as.precintcon.deciles.r
#' @include deciles.r
NULL

#' @name as.deciles
#' @author Lucas Venezian Povoa \email{lucasvenez@@gmail.com} 
#' @aliases as.precintcon.deciles as.deciles deciles precintcon.deciles.analysis
#' @title Deciles of a precipitation serie 
#' @description Groups the monthly precipitation into decis, 
#' i.e., it splits a precipitation serie into ten equal parts in crescent 
#' order, from the lower to the highest precipitation. 
#' @usage as.precintcon.deciles(object) 
#' @param object a daily or monthly precipitation serie.
#' @return A data.frame (precintcon.deciles) containing the following variables:
#' \itemize{
#' \item \code{D1} corresponds to the precipitation values not exceeding 10\% 
#' of the lowest values.
#' \item \code{D2} corresponds to the precipitation values not exceeding 20% 
#' of the lowest values.
#' \item \code{D3} corresponds to the precipitation values not exceeding 30% 
#' of the lowest values.
#' \item \code{D4} corresponds to the precipitation values not exceeding 40\% 
#' of the lowest values.
#' \item \code{D5} is equals to the median that corresponds to the precipitation 
#' values not exceeding 50\% of the lowest values.
#' \item \code{D6} corresponds to the precipitation values not exceeding 60\% 
#' of the lowest values.
#' \item \code{D7} corresponds to the precipitation values not exceeding 70% 
#' of the lowest values.
#' \item \code{D8} corresponds to the precipitation values not exceeding 80% 
#' of the lowest values.
#' \item \code{D9} corresponds to the precipitation values not exceeding 90% 
#' of the lowest values.
#' \item \code{D10} corresponds to the precipitation values not exceeding 100% 
#' of the lowest values.
#' } 
#' @seealso 
#' \code{\link{read.data}} 
#' @keywords precintcon deciles
#' @export 
as.deciles <- function(object) {
   return(as.precintcon.deciles(object))  
}