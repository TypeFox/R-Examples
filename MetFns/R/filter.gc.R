filter.gc<-function(data,long.low=0,long.up=180, ew=c("E","W"), lat.low=0,lat.up=90, ns=c("N","S"))
{
   if(!is.data.frame(data) || !is.numeric(c(lat.low,lat.up,long.low,long.up)) || !is.character(c(ns,ew)) || 
       (any(c(lat.low,long.low)<0) || lat.up>90 || long.up>180) ||  long.low>long.up || lat.low>lat.up || 
       !ns%in%c("N","S") || !ew%in%c("E","W") ) 
      stop("invalid input parameter(s) specification: check data/long.low/long.up/ew/lat.low/lat.up/ns")
  
   data[data$long>=long.low & data$long<=long.up & data$EW==ew & data$lat>=lat.low & data$lat<=lat.up & data$NS==ns,]
}
