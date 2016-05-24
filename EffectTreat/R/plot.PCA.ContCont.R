plot.PCA.ContCont <- plot.Multivar.PCA.ContCont <- function(x, Xlab.PCA, Main.PCA, Type="Percent", 
                              Labels=FALSE, PCA=TRUE, 
                              Good.Pretreat=FALSE, EffectT0T1=FALSE, Main.Good.Pretreat, 
                              Par=par(oma=c(0, 0, 0, 0), mar=c(5.1, 4.1, 4.1, 2.1)), col, ...){   

  Object <- x 
  
  if (methods::is(Object, "Multivar.PCA.ContCont")==TRUE & Good.Pretreat==TRUE){
    stop("A plot of delta (argument Good.Pretreat=TRUE can only be obtained for fitted objects of class PCA.ContCont.")
  }
  
  
  # voor uni
  if (methods::is(Object, "PCA.ContCont")==TRUE){
  
  if (missing(Xlab.PCA)) {Xlab.PCA <- expression(rho[psi])}
  if (missing(Main.PCA)) {Main.PCA="PCA"} 
  if (missing(col)) {col <- c(8)}
  
  if (PCA==TRUE){
    
    dev.new() 
    par=Par  
    if (Type=="Freq"){
      h <- hist(Object$PCA, ...)
      h$density <- h$counts/sum(h$counts)
      cumulMidPoint <- ecdf(x=Object$PCA)(h$mids)
      labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
      
      if (Labels==FALSE){
        plot(h,freq=T, xlab=Xlab.PCA, ylab="Frequency", col=col, main=Main.PCA, ...)
      }
      if (Labels==TRUE){
        plot(h,freq=T, xlab=Xlab.PCA, ylab="Frequency", col=col, main=Main.PCA, labels=labs, ...)
      }
    }
    
    if (Type=="Percent"){
      h <- hist(Object$PCA, ...)
      h$density <- h$counts/sum(h$counts)
      cumulMidPoint <- ecdf(x=Object$PCA)(h$mids)
      labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
      
      if (Labels==FALSE){
        plot(h,freq=F, xlab=Xlab.PCA, ylab="Percentage", col=col, main=Main.PCA, ...)
      }
      if (Labels==TRUE){
        plot(h,freq=F, xlab=Xlab.PCA, ylab="Percentage", col=col, main=Main.PCA, labels=labs, ...)
      }
    }
    
    if (Type=="CumPerc"){
      h <- hist(Object$PCA, breaks=length(Object$PCA), ...)
      h$density <- h$counts/sum(h$counts)
      cumulative <- cumsum(h$density)
      plot(x=h$mids, y=cumulative, xlab=Xlab.PCA, ylab="Cumulative percentage", col=0, main=Main.PCA, ...)
      lines(x=h$mids, y=cumulative)
    }
  }
  
  
  if (Good.Pretreat==TRUE){
    
    if (missing(Main.Good.Pretreat)) {Main.Good.Pretreat = " "}  
    par=Par
    if (Type=="Freq"){
      h <- hist(Object$GoodSurr$delta, ...)
      h$density <- h$counts/sum(h$counts)
      cumulMidPoint <- ecdf(x=Object$GoodSurr$delta)(h$mids)
      labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
      
      if (Labels==FALSE){
        plot(h,freq=T, xlab=expression(delta), ylab="Frequency", main=Main.Good.Pretreat, col=col, ...)
      }
      if (Labels==TRUE){
        plot(h,freq=T, xlab=expression(delta), ylab="Frequency", col=col, labels=labs, 
             main=Main.Good.Pretreat, ...)
      }
    }
    
    if (Type=="Percent"){
      h <- hist(Object$GoodSurr$delta, ...)
      h$density <- h$counts/sum(h$counts)
      cumulMidPoint <- ecdf(x=Object$GoodSurr$delta)(h$mids)
      labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
      
      if (Labels==FALSE){
        plot(h,freq=F, xlab=expression(delta), ylab="Percentage", col=col, main=Main.Good.Pretreat, ...)
      }
      if (Labels==TRUE){
        plot(h,freq=F, xlab=expression(delta), ylab="Percentage", col=col, labels=labs, 
             main=Main.Good.Pretreat, ...)
      }
    }
    
    if (Type=="CumPerc"){
      h <- hist(Object$GoodSurr$delta, breaks=length(Object$GoodSurr$delta), ...)
      h$density <- h$counts/sum(h$counts)
      cumulative <- cumsum(h$density)
      plot(x=h$mids, y=cumulative, xlab=expression(delta), ylab="Cumulative percentage", 
           col=0, main=Main.Good.Pretreat, ...)
      lines(x=h$mids, y=cumulative)
    }
    
    
  }    
  
  if (EffectT0T1==TRUE){
  
    if (methods::is(Object, "PCA.ContCont")==TRUE){
    if (length(unique(Object$Pos.Def$T0S))>1){
    cat("WARNING: This plot shows PCA as a function of rho_T0T1. Rho_T0S is not constant here, so this plot may be misleading")  
    }
    if (length(unique(Object$Pos.Def$T1S))>1){
      cat("WARNING: This plot shows PCA as a function of rho_T0T1. Rho_T1S is not constant here, so this plot may be misleading")  
    }
    
    plot(Object$Pos.Def$T0T1, Object$PCA, xlab=expression(rho[T0T1]), 
         ylab=expression(rho[psi]), col=0, ...)
    lines(Object$Pos.Def$T0T1, Object$PCA)
    }
  }
    
  } # einde univar

    
    
  # voor multi
    if (methods::is(Object, "Multivar.PCA.ContCont")==TRUE){  
      if (missing(Xlab.PCA)) {Xlab.PCA <- expression(R[psi]^2)}
      if (missing(Main.PCA)) {Main.PCA="PCA"} 
      if (missing(col)) {col <- c(8)}
      
      if (PCA==TRUE){
        
        dev.new() 
        par=Par  
        if (Type=="Freq"){
          h <- hist(Object$PCA, ...)
          h$density <- h$counts/sum(h$counts)
          cumulMidPoint <- ecdf(x=Object$PCA)(h$mids)
          labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
          
          if (Labels==FALSE){
            plot(h,freq=T, xlab=Xlab.PCA, ylab="Frequency", col=col, main=Main.PCA, ...)
          }
          if (Labels==TRUE){
            plot(h,freq=T, xlab=Xlab.PCA, ylab="Frequency", col=col, main=Main.PCA, labels=labs, ...)
          }
        }
        
        if (Type=="Percent"){
          h <- hist(Object$PCA, ...)
          h$density <- h$counts/sum(h$counts)
          cumulMidPoint <- ecdf(x=Object$PCA)(h$mids)
          labs <- paste(round((1-cumulMidPoint), digits=4)*100, "%", sep="")
          
          if (Labels==FALSE){
            plot(h,freq=F, xlab=Xlab.PCA, ylab="Percentage", col=col, main=Main.PCA, ...)
          }
          if (Labels==TRUE){
            plot(h,freq=F, xlab=Xlab.PCA, ylab="Percentage", col=col, main=Main.PCA, labels=labs, ...)
          }
        }
        
        if (Type=="CumPerc"){
          h <- hist(Object$PCA, breaks=length(Object$PCA), ...)
          h$density <- h$counts/sum(h$counts)
          cumulative <- cumsum(h$density)
          plot(x=h$mids, y=cumulative, xlab=Xlab.PCA, ylab="Cumulative percentage", col=0, main=Main.PCA, ...)
          lines(x=h$mids, y=cumulative)
        }
      }
      
      
      if (EffectT0T1==TRUE){
          plot(Object$T0T1, Object$PCA, xlab=expression(rho[T0T1]), 
               ylab="PCA", col=0, ...)
          lines(Object$T0T1, Object$PCA)
      }
    }
  
}