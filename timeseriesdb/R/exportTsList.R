#' Export A List of Time Series to CSV
#' 
#' Export a List of time series to semi-colon separated csv file. Typically multiple time series are exported as long format (melted format) files. 
#' 
#' @param tl list of time series
#' @param fname character file name. If set to NULL a standard file name chunk + Sys.Date is used.
#' @param cast logical. Should the resulting data.frame be cast to wide format? Defaults to TRUE
#' @importFrom reshape2 dcast 
#' @export
exportTsList <- function(tl,fname = NULL,cast = T){
  tl <- lapply(tl,as.xts)
  out_list <- lapply(names(tl),function(x){
    dframe <- data.frame(time = time(tl[[x]]),
                         value = tl[[x]],row.names = NULL)
    dframe$series <- x
    dframe
  })
  
  tsdf <- do.call("rbind",out_list)
#   if(is.null(fname)) fname <- "timeseriesdb_export_"
#   fname <- paste0(fname,gsub("-","_",Sys.Date()),".csv")
  #write.csv2(tsdf,file = fname,row.names = F)
#  reshape2::dcast(tsdf)
  if(cast){
    tsdf <- reshape2::dcast(tsdf,time ~ series)
  }
  
  write.csv2(tsdf,file = fname,row.names = F)              
}
