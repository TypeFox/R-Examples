filter.totcor<-function(data,shw,Ralpha=NULL,Delta=NULL,r,C=5)
{
 if(!is.numeric(c(r,C)) || (C<1)) 
      stop("invalid input parameter(s) specification: check r/C")

 elev.data<-sinh(data,shw, Ralpha, Delta)
 m<-ncol(elev.data)
 elev.data[elev.data$F*r^(6.5-elev.data$lmg)/elev.data$sine.h<C,-m]
}
