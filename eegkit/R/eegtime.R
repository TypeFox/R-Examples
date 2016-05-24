eegtime <- 
  function(time,voltage,flipvoltage=TRUE,vlty=1,vlwd=2,vcol="blue",
           voltageSE=NULL,slty=NA,slwd=1,scol="cyan",salpha=0.65,conflevel=0.95,
           plotzero=TRUE,zlty=1,zlwd=0.5,zcol="black",xlim=NULL,ylim=NULL,
           xlab=NULL,ylab=NULL,nxtick=6,nytick=6,add=FALSE,...){
    ###### Plots Single-Channel EEG Time Course
    ###### Nathaniel E. Helwig (helwig@umn.edu)
    ###### Last modified: February 16, 2015
    
    ### initial checks
    vlen <- length(voltage)
    if(length(time)!=vlen){stop("Inputs 'time' and 'voltage' must be same length.")}
    if(!is.null(xlim[1])){
      if(length(xlim)!=2L){stop("Input 'xlim' must be two element vector.")}
      if(xlim[1]>min(time)){stop("First input of 'xlim' is larger than minimum 'time' input.")}
      if(xlim[2]<max(time)){stop("Second input of 'xlim' is smaller than maximum 'time' input.")}
    } else {
      xlim <- range(time)
    }
    if(!is.null(ylim[1])){
      if(length(ylim)!=2L){stop("Input 'ylim' must be two element vector.")}
      if(ylim[1]>min(voltage)){stop("First input of 'ylim' is larger than minimum 'voltage' input.")}
      if(ylim[2]<max(voltage)){stop("Second input of 'ylim' is smaller than maximum 'voltage' input.")}
    } else {
      ylim <- range(voltage)
    }
    if(is.null(xlab[1])){xlab <- "Time After Stimulus (ms)"}
    if(is.null(ylab[1])){ylab <- expression("Voltage ("*mu*"V)")}
    nxtick <- as.integer(nxtick[1])
    nytick <- as.integer(nytick[1])
    if(any(c(nxtick,nytick)<1)){stop("Inputs 'nxtick' and  'nytick' must be positive integers.")}
    
    ### check for standard errors
    if(is.null(voltageSE)==FALSE){
      if(length(voltageSE)!=vlen){stop("Inputs 'voltage' and 'voltageSE' must be same length.")}
      if(any(voltageSE<0)){stop("Input 'voltageSE' must contain nonnegative standard errors.")}
      if(conflevel<=0 || conflevel>=1){stop("Input 'conflevel' must statisfy 0<conflevel<1.")}
    }
    
    ### make x and y tick marks
    xticks <- pretty(seq(xlim[1],xlim[2],l=nxtick))
    yticks <- pretty(seq(ylim[1],ylim[2],l=nytick))
    
    ### possibly flip voltage and yticks
    if(flipvoltage){
      voltage <- (-1)*voltage
      ylim <- (-1)*rev(ylim)
      yticks <- (-1)*rev(yticks)
      ytlabs <- abs(yticks)
      for(ii in 1:length(yticks)){
        if(yticks[ii]<0){ytlabs[ii] <- paste("+",ytlabs[ii],sep="")}
        if(yticks[ii]>0){ytlabs[ii] <- paste("-",ytlabs[ii],sep="")}
      }
    } else {
      ytlabs <- yticks
    }
    
    ### plot voltage
    if(add){
      lines(time,voltage,lwd=vlwd,lty=vlty,col=vcol,...)
    } else {
      plot(time,voltage,type="l",lwd=vlwd,lty=vlty,col=vcol,
           xlab=xlab,ylab=ylab,axes=FALSE,xlim=xlim,ylim=ylim,...)
      axis(1,at=xticks,...)
      axis(2,at=yticks,labels=ytlabs,...)
    }
    
    ### plot voltageSE and zero line
    if(!is.null(voltageSE)){
      cval <- qnorm(1-(1-conflevel)/2)
      if(is.na(slty)){
        myrgb <- col2rgb(scol)/255
        polycoords <- rbind(cbind(time,voltage+cval*voltageSE),
                            cbind(rev(time),rev(voltage-cval*voltageSE)))
        polygon(polycoords[,1],polycoords[,2],col=rgb(myrgb[1],myrgb[2],myrgb[3],salpha),border=NA)
        lines(time,voltage,lwd=vlwd,lty=vlty,col=vcol)
      } else {
        lines(time,voltage-cval*voltageSE,lwd=slwd,lty=slty,col=scol)
        lines(time,voltage+cval*voltageSE,lwd=slwd,lty=slty,col=scol)
      }
    }
    if(plotzero){lines(xlim,rep(0,2),lwd=zlwd,lty=zlty,col=zcol)}
    
}