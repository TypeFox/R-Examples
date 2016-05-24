#' Import A Long Format CSV File Into a List of Time Series
#' 
#' Generates a List of time series of class tslist from a flat file that contains multiple time series. 
#' 
#' @param f character file name
#' @param base_ts logical should the export be casted to R basic time series or left as xts? Defaults to TRUE.  
#' @importFrom zoo as.yearmon
#' @importFrom zoo zoo 
#' @importFrom xts as.xts
#' @export 
importTsList <- function(f,base_ts = T){
  dframe <- read.csv2(f,header = T,sep = ";",
                      stringsAsFactors = F)
  l <- split(dframe,factor(dframe$series))
  l <- lapply(l,function(x){
    xts::as.xts(x$value,zoo::as.yearmon(x$time))
  })
  
  if(base_ts){
    l <- lapply(l,function(x) as.ts(zoo::zoo(x)))
  }
  
  class(l) <- append(class(l),"tslist")
  l
}
