filter.h<-function(data,shw,Ralpha=NULL,Delta=NULL,h.low=10,h.up=90)
{
 if(!is.numeric(c(h.low,h.up)) || (h.low<10 || h.up>90) || h.low>h.up) 
      stop("invalid input parameter(s) specification: check h.low/h.up")

 elev.data<-sinh(data,shw, Ralpha, Delta)
 m<-ncol(elev.data)
 h<-asin(elev.data$sine.h)*180/pi
 elev.data[h>=h.low & h<=h.up,-m]
}
