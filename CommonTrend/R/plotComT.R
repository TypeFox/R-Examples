plotComT <-function(ComT, i, 
                   x.axis=NA, approx.ticks=7,
                   legend=c("Original Data","Common Trend"), 
                   main="", ylab="", xlab="") 
  
{  if(!class(ComT)=="ComT") stop ("Class must be 'ComT' ")

   ComTD=ComT$common.trend[i,]
   Origin=ComT$data.used[i,]
   max=length(ComTD) 
   
   plot (Origin,type="l",ylab=ylab,xlab=xlab, axes  =  FALSE)
   
   if (length(x.axis) == (max+ComT$lag.chosen) ) { 
                 xlable=paste (x.axis[(length(x.axis)-max+1):length(x.axis)])
                 ticks=seq(1,max, by=round(max/approx.ticks))
                 axis(1,at=ticks,labels= xlable[ticks], cex.axis=0.7) 
   }else if (is.na(x.axis[1])==TRUE) {
                 axis(1)
   }else if (!length(x.axis) == (max+ComT$lag.chosen) ){ 
                         stop("Length of x.axis must be the same as length of original data" ) 
          } 
           
   axis(2)
   title(main=main)
   lines (ComTD+mean(ComT$stationary[i,]),lwd=2,col="BLUE") 
   
   position=max(Origin)
   legend(0,position,legend,lwd=1:2,col=c("Black","BLUE"),bty ="n")
   box()
}

