#' @name pci.supraseasonal
#' @aliases pci.supraseasonal
#' @author Lucas Venezian Povoa \email{lucasvenez@@gmail.com}
#' 
#' @title Supraseasonal Precipitation Concentration Index
#' @description Calculates the supraseasonal precipitation concentration index.
#' @param object is a daily or monthly precipitation serie.
#' @param hemisthere is the hemisthere, "n" for northern and "s" for south, of the precipitation serie.
#' @return A data.frame containing the following variables:
#' \itemize{
#' \item \code{year} is the year;
#' \item \code{season} is the meteorological supraseason, wet or dry; and
#' \item \code{pci.season} is the seasonal perceptation concentration index.
#' }
#' @examples
#' ##
#' # Loading the daily precipitation serie
#' data(daily)
#' 
#' ##
#' # Calculating the supraseasonal precipitation concentration index
#' pci.supraseasonal(daily, hemisthere = "s")
#' @references
#' M. de Luis, J. C. Gonz\'alez-Hidalgo, M. Brunetti,  L. A. Longares (2011). 
#' Precipitation concentration changes in Spain 1946-2005. Natural Hazards and Earth System Science, 
#' 5:11, pp. 1259--1265
#' @export
pci.supraseasonal <- function(object, hemisthere = c("n", "s")) {
  
  object <- as.precintcon.monthly(object)
  
  result <- data.frame()
  
  station <- c("dry", "wet")
  
  start <- which(object$month == 4)[1]
  
  for(i in seq(start, nrow(object), by = 6)) {
    
    if (nrow(object) - i < 5) break;
    
    pci <- (sum(object[i:(i+5),3]**2)/sum(object[i:(i+5),3])**2)*50

    pci <- ifelse(is.nan(pci), 0, pci)
    
    result <- rbind(result, 
              data.frame(
                year              = object[i,1], 
                season            = station[abs((i - start) %/% 6 + ifelse(hemisthere == "n", 0, 1)) %% 2 + 1], 
                pci.supraseasonal = pci))
  }
  
  class(result) <- c("data.frame", "precintcon.seasonal")
  
  return(result)
}