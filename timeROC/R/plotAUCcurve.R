plotAUCcurve <- function(object,FP=2,add=FALSE,conf.int=FALSE,conf.band=FALSE,col="black"){
  ## browser()
  # {{{ class "ipcwsurvivalROC"
  if(class(object)=="ipcwsurvivalROC"){
    AUC <-  object$AUC[!is.na(object$AUC)]
    times <- object$times[!is.na(object$AUC)]
    if (add==FALSE){
      plot(times,AUC,xlab="time t",ylab="AUC(t)",lwd=2,ylim=c(0.45,1),type="l",col=col)
      abline(h=0.5,lty=2)
    }else{
      lines(times,AUC,lwd=2,type="l",col=col)
    }    
    if(object$iid==TRUE & (conf.int==TRUE | conf.band==TRUE)){
      result.confint <- confint(object=object,level=0.95,n.sim=3000)
    }
    if(conf.int==TRUE){
      if(object$iid==TRUE){
        matlines(times,result.confint$CI_AUC/100,lty=2,lwd=2,col=col)
      }else{
        cat("Confidence intervals cannot be computed because you have chosen iid=FALSE (default) in the input of timeROC(). \n")
      }
    }
    if(conf.band==TRUE){
      if(object$iid==TRUE){
        matlines(times,result.confint$CB_AUC/100,lty=3,lwd=2,col=col)
      }else{
        cat("Confidence band cannot be computed because you have chosen iid=FALSE (default) in the input of timeROC(). \n")
      }
    }
  }
  # }}}
  # {{{ class "ipcwcompetingrisksROC"
  if(class(object)=="ipcwcompetingrisksROC"){
    if (FP==2){
      AUC <-  object$AUC_2[!is.na(object$AUC_2)]
      times <- object$times[!is.na(object$AUC_2)]
    }
    if (FP==1){
      AUC <-  object$AUC_1[!is.na(object$AUC_1)]
      times <- object$times[!is.na(object$AUC_1)]
    }
    if (FP!=1 & FP!=2){
      stop("FP indicates the type of specificity. It must be 1 or 2 \n")
    }
    if (add==FALSE){
      plot(times,AUC,xlab="time t",ylab="AUC(t)",lwd=2,ylim=c(0.45,1),type="l",col=col)
      abline(h=0.5,lty=2)
    }else{
      lines(times,AUC,lwd=2,type="l",col=col)
    }    
    if(object$iid==TRUE & (conf.int==TRUE | conf.band==TRUE)){
      result.confint <- confint(object=object,level=0.95,n.sim=3000)
    }
    if(conf.int==TRUE){
      if(object$iid==TRUE){
        if (FP==2){
          matlines(times,result.confint$CI_AUC_2/100,lty=2,lwd=2,col=col)
        }
        if (FP==1){
          matlines(times,result.confint$CI_AUC_1/100,lty=2,lwd=2,col=col)
        }
      }else{
        cat("Confidence intervals cannot be computed because you have chosen iid=FALSE (default) in the input of timeROC(). \n")
      }
    }
    if(conf.band==TRUE){
      if(object$iid==TRUE){
        if (FP==2){
          matlines(times,result.confint$CB_AUC_2/100,lty=3,lwd=2,col=col)
        }
        if (FP==1){
          matlines(times,result.confint$CB_AUC_1/100,lty=3,lwd=2,col=col)
        }
      }else{
        cat("Confidence band cannot be computed because you have chosen iid=FALSE (default) in the input of timeROC(). \n")
      }
    }
  }
  # }}}
  # {{{ error if not class ipcwcompetingrisksROC nor ipcwsurvivalROC
  if(class(object)!="ipcwcompetingrisksROC" & class(object)!="ipcwsurvivalROC"){
    stop("Function written for object of class ipcwsurvivalROC or ipcwcompetingrisksROC\n")
  }
  # }}}
}
