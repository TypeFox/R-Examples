plot.ICA.BinCont <- function(x, Xlab, Main=" ", Type="Density", Labels=TRUE, ...){

  if (missing(Xlab)) {Xlab <- expression(R[H]^2)}
  
  Object <- x 
  
  dev.new()
  
    if (Type=="Density"){    
    plot(density(x$R2_H, na.rm = T), xlab=Xlab, ylab="Density", main=Main, lwd=2)
    }
    
    
    if (Type=="Freq"){
      h <- hist(Object$R2_H, plot = FALSE, ...)
      h$density <- h$counts/sum(h$counts)
      cumulMidPoint <- ecdf(x=Object$R2_H)(h$mids)
      labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
      
      if (Labels==FALSE){
        plot(h,freq=T, xlab=Xlab, ylab="Frequency", main=Main)
      }
      if (Labels==TRUE){
        plot(h,freq=T, xlab=Xlab, ylab="Frequency", main=Main, labels=labs)
      }
    }
    
    
    if (Type=="Percent"){
      
      h <- hist(Object$R2_H, plot = FALSE, ...)
      h$density <- h$counts/sum(h$counts)
      cumulMidPoint <- ecdf(x=Object$R2_H)(h$mids)
      labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
      
      if (Labels==FALSE){
        plot(h, freq=F, xlab=Xlab, ylab="Percentage", main=Main)
      }
      if (Labels==TRUE){
        plot(h, freq=F, xlab=Xlab, ylab="Percentage", main=Main, labels=labs)
      }
      }
      
  
    if (Type=="CumPerc"){
      
      h <- hist(Object$R2_H, breaks=length(Object$R2_H), plot = FALSE, ...)
      h$density <- h$counts/sum(h$counts)
      cumulative <- cumsum(h$density)
      plot(x=h$mids, y=cumulative, xlab=Xlab, ylab="Cumulative percentage", col=0, main=Main)
      lines(x=h$mids, y=cumulative)
    }
  }

 