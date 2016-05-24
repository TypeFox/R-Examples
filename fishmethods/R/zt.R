zt<-function(x=NULL,ctype=1){
    if(is.null(x)) stop ("time series is missing.")
    if(!is.numeric(x)) stop ("time series is not numeric.")
    if(ctype==1) return((x-mean(x,na.rm=T))/sd(x,na.rm=T))
    if(ctype==2) return(x-mean(x,na.rm=T))
   }
