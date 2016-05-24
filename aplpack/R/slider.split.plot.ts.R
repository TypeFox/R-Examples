slider.split.plot.ts<-function(x,type="l",...)
{  
 x.name<-deparse(substitute(x))
 if(missing(x)||length(x)<2) return("Error: x must be a vector")
 y<-x; if(is.ts(x)) { x<-time(x) } else { x<-seq(x) }

 args<-list(...)
 n<-length(x); xmin<-min(x); xmax<-max(x)
 xdelta<-xmax-xmin
 slider(obj.name="summary.type",obj.value="linear")
 refresh<-function(...){
 # initialization
   summary.type<-slider(obj.name="summary.type")
   width<-slider(no=1)
   limit<-slider(no=2)
   n.sec<-1 
   limit<-limit-width*ceiling((limit-xmin)/width)
 # plot: # plot(x,y,type=type,bty="n",xlab="",ylab="")
   do.call("plot",c(alist(x,y,type=type),args))
   limit<-limit-width-width/n.sec; j<-0
 #   abline(v=limits,lwd=0.5,lty=3)    
   while(limit<xmax){ j<-j+1
     limit<-limit+width/n.sec; limit2<-limit+width
     ind<-limit<=x & x<=limit2
     xx<-x[ind]; yy<-y[ind]; if(length(xx)<2) next
     abline(v=limit,lwd=0.5,lty=3,col=1)    
     if(summary.type=="linear"){
       coef<-lm(yy~xx)$coef
       segments(limit, coef[1]+coef[2]*limit,
                limit2,coef[1]+coef[2]*limit2,col=j) 
     }
     if(summary.type=="five.number"){
       five<-fivenum(yy)
       xx<-0.5*(limit+limit2)
       points(xx,five[3],pch=19,cex=1,col="red")
       segments(xx,five[1],xx,five[2],lwd=3,col="red")
       segments(xx,five[4],xx,five[5],lwd=3,col="red")
     }
   }
   abline(v=slider(no=2),lwd=0.5,lty=1,col=1)    
 }    
 f1<-function(...){
   slider(obj.name="summary.type",obj.value="linear")
   refresh()
 }
 f2<-function(...){
   slider(obj.name="summary.type",obj.value="five.number")
   refresh()
 }
 slider(refresh,c("width of window","limit"),
      c(xdelta/length(x)*3,xmin),c(xdelta,xmax),
      c(xdelta/1000,xdelta/1000),c(xdelta/4,xmin),
      but.functions=c(f1,f2),
      but.names=c("linear model","fivenum summary")
 )
 refresh() 
 cat("select window and summary type and look at time series!\n")
}

