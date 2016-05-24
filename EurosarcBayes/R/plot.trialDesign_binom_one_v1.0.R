





.plot.trialDesign_binom_one=function(x,xaxs="i",yaxs="i",xlim=NULL,ylim=NULL,lwd=2,cols=NULL,lty=NULL,col.s=3,col.f=2,main=NULL,...){


  # single design provided
  n=as.numeric(x@reviews)
  options(warn=-1)
  success=as.numeric(x@success)
  failure=as.numeric(x@failure)
  options(warn=0)

  if(is.null(xlim)){
    xlim=c(0,max(n)+2)
  }
  if(is.null(ylim)){
    ylim=c(0,max(success,na.rm = TRUE)+2)
  }


  plot(c(0,n),c(0,success),type="n",xlim=xlim,ylim=ylim,xlab="Number of patients recruited",ylab="Number of successes",xaxs=xaxs,yaxs=yaxs,...)

  if(length(success)>0){
    points(n,success,col=col.s,pch=16)
    if(length(success)<10){
      for(i in 1:length(success)){
        if(!is.na(success[i])){
          segments(n[i],success[i],n[i],n[i],col=col.s,lwd=lwd)
        }
      }
    }
  }

  if(length(failure)>0){
    points(n,failure,col=col.f,pch=16)
    if(length(failure)<10){
      for(i in 1:length(failure)){
        if(!is.na(failure[i])){
          segments(n[i],failure[i],n[i],0,col=col.f,lwd=lwd)
        }
      }
    }
  }
}

setMethod(f="plot",signature=c(x="trialDesign_binom_one"),definition=.plot.trialDesign_binom_one)


# ended
