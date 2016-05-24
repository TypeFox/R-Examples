date2winter<-function(x,first=10,last=4){
  if(!first%in%1:12) stop("'first' must be a month number")
  if(!last%in%1:12) stop("'last' must be a month number")
  if (inherits(x,"POSIXct")) x<-as.POSIXlt(x) 
  if (!inherits(x,"POSIXlt")) stop("x must be of class POSIXct or POSIXlt")
  res<-ifelse((x$mon+1)>last & (x$mon+1)<first,"Excluded",ifelse((x$mon+1)>=first,paste(x$year+1900,"-",x$year+1900+1,sep=""),paste(x$year+1900-1,"-",x$year+1900,sep="")))
  attributes(res)$span<-c(first=first,last=last)
  res
}


