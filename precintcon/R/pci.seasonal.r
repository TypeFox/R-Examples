#' @name pci.seasonal
#' @aliases pci.seasonal
#' @author Lucas Venezian Povoa \email{lucasvenez@@gmail.com}
#' @description Calculates the Precipitation Concentration Index (PCI) in a seasonal granularity 
#' on a daily or monthly precipitation serie.
#' @title Seasonal Precipitation Concentration Index
#' @param object is a daily or monthly precipitation serie
#' @param hemisthere is the hemisthere, "n" for northern and "s" for south, of the precipitation serie
#' @return A data.frame containing the following variables:
#' \itemize{
#' \item \code{year} is the year;
#' \item \code{season} is the meteorological season; and
#' \item \code{pci.seasonal} is the seasonal perceptation concentration index.
#' }
#' @examples 
#' ##
#' # Loading the daily precipitation serie
#' data(daily)
#' 
#' ##
#' # Calculating the seasonal perceptation concentration index
#' pci.seasonal(daily, hemisthere = "s")
#' @export
pci.seasonal <- function(object, hemisthere) {
  
  object <- as.precintcon.monthly(object)
  
  result <- data.frame()
  
  station <- c("spring", "summer", "autumn", "winter")
  
  start <- which(object$month == 3)[1]
  
  for(i in seq(start, nrow(object), by = 3)) {
    
    if (nrow(object) - i < 2) break
    
    pci <- (sum(object[i:(i+2),3]**2)/sum(object[i:(i+2),3])**2)*25
    
    pci <- ifelse(is.nan(pci), 0, pci)
    
    result <- rbind(result, 
              data.frame(
                year         = object[i,1], 
                season       = station[abs((i - start) %/% 3 + ifelse(hemisthere == "n", 0, 2)) %% 4 + 1], 
                pci.seasonal = pci))
  }
  
  class(result) <- c("data.frame", "precintcon.seasonal")
  
  return(result)
}