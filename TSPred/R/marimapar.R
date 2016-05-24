marimapar <-
function(TimeSeries, na.action=na.omit, xreg=NULL){
  return (lapply(TimeSeries,arimapar,na.action=na.action, xreg=xreg))
}
