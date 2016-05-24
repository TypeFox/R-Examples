lag2.plot <-
function(series1,series2,max.lag=0,corr=TRUE,smooth=TRUE){ 
   name1=paste(deparse(substitute(series1)),"(t-",sep="")
   name2=paste(deparse(substitute(series2)),"(t)",sep="")
   series1=as.ts(series1)
   series2=as.ts(series2)
   max.lag=as.integer(max.lag)
   m1=max.lag+1
   prow=ceiling(sqrt(m1))
   pcol=ceiling(m1/prow)
   a=stats::ccf(series1,series2,max.lag,plot=FALSE)$acf
   old.par <- par(no.readonly = TRUE)
   par(mfrow=c(prow,pcol), mar=c(2.5, 4, 2.5, 1), cex.main=1.1, font.main=1)
   for(h in 0:max.lag){                   
   plot(lag(series1,-h), series2, xy.labels=FALSE, main=paste(name1,h,")",sep=""), ylab=name2, xlab="") 
    if (smooth==TRUE) 
    lines(stats::lowess(ts.intersect(lag(series1,-h),series2)[,1],
                 ts.intersect(lag(series1,-h),series2)[,2]), col="red")
    if (corr==TRUE)
    legend("topright", legend=round(a[m1-h], digits=2), text.col ="blue", bg="white", x.intersp=0)             
   on.exit(par(old.par))
   }
}

