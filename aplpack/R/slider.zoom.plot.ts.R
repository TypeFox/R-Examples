slider.zoom.plot.ts<-function(x,n.windows=1,...)
{ 
 x.name<-deparse(substitute(x))
 if(missing(x)||length(x)<2) return("Error: x must be a vector")
 y<-x; if(is.ts(x)) { x<-time(x) } else { x<-seq(x) }

 args<-list(...)
 tmin<-1; tmax<-length(x)
 refresh<-function(...){
 # initialization
   width <- slider(no=1); tstart1 <- slider(no=2); tend1 <- width+tstart1
   if(n.windows>1){
     tstart2 <- slider(no=3); tend2 <- width+tstart2
   }
 # plot
   par(mfrow=c(2+(n.windows>1),1),mai=c(.5,0.5,.1,0))
   # plot(x,y,type="l",bty="n",xlab="",ylab="")
   do.call("plot",c(alist(x,y,type="l",bty="n"),args))
   abline(v=c(x[tstart1],x[tend1]),col="red")
   lines(x[tstart1:tend1],y[tstart1:tend1],col="red",lty=2)
   if(n.windows>1){
     abline(v=c(x[tstart2],x[tend2]),col="blue")
     lines(x[tstart2:tend2],y[tstart2:tend2],col="blue",lty=3)
   }
   usr<-par()$usr
   ind<-tstart1:tend1
   plot(x[ind],y[ind],type="b",col="red",bty="n", # ylim=usr[3:4],
        xlim=c(x[tstart1],x[tstart1]+width*diff(x[1:2])))
   if(n.windows>1){
   ind<-tstart2:tend2
     plot(x[ind],y[ind],type="b",col="blue",bty="n", # ylim=usr[3:4],
        xlim=c(x[tstart2],x[tstart2]+width*diff(x[1:2])))
   }
   par(mfrow=c(1,1))
 }
 if(n.windows<2){
   slider(refresh,c("width of window","start of window"),
          c(1,1),c(tmax,tmax),c(1,1),c(ceiling(tmax/4),1))
 }else{
   slider(refresh,
      c("width of window","start window 1","start window 2"),
      c(1,1,1),c(tmax,tmax,tmax),c(1,1,1),c(ceiling(tmax/4),1,ceiling(tmax/2)))
 }
 refresh() 
 cat("select window and look at time series!\n")
}

