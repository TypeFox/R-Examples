 gap<-function(x=NULL){
         if(is.null(x)) stop ("x is missing")
         if(!is.numeric(x)) stop ("data (x) are not numeric.")
        n<-length(x)
        x<-sort(x)
        z<-NULL
        for(i in 1:c(n-1)) z[n-i+1]<-sqrt((i*(n-i)*(x[n-i+1]-x[n-i])))
        zt<-mean(z[-c(1)],trim=0.25,na.rm=T)
        z/zt
  }
 