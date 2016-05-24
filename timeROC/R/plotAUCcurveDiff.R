plotAUCcurveDiff <- function(object1,object2,FP=2,add=FALSE,conf.int=FALSE,conf.band=FALSE,col="black",ylim=c(-0.5,0.5)){
 ## browser()
  # {{{ check class of object1 and 2 are equal
  if (class(object1)!=class(object2)){
    stop("object1 and object2 must be of the same class \n")
  }
  if(!identical(object1$times,object2$times)){
    stop(paste("The two objects you want to compare should have been computed for the same vector of times\n",sep=""))
  }
  if(!object1$n==object2$n){
    stop(paste("The two objects you want to compare were not fitted with the same subjects.\n
               This function does not deal with such a case.",sep=""))
  } 
  # }}}
  # {{{ class "ipcwsurvivalROC"
  if(class(object1)=="ipcwsurvivalROC"){
    AUC <-  object1$AUC[!(is.na(object1$AUC)|is.na(object2$AUC))]- object2$AUC[!(is.na(object1$AUC)|is.na(object2$AUC))]
    times <- object1$times[!(is.na(object1$AUC)|is.na(object2$AUC))]
    if (add==FALSE){
      plot(times,AUC,xlab="time t",ylab=expression(paste(Delta,"AUC(t)")),lwd=2,ylim=ylim,type="l",col=col)
      abline(h=0,lty=2)
    }else{
      lines(times,AUC,lwd=2,type="l",col=col)
    }
    if(object1$iid==TRUE & object2$iid==TRUE & (conf.int==TRUE | conf.band==TRUE)){
      mat.iid <- object1$inference$mat_iid_rep_1-object2$inference$mat_iid_rep_1
      se<-apply(mat.iid,2,sd)/sqrt(object1$n)
    }
    if(conf.int==TRUE){
      if(object1$iid==TRUE & object2$iid==TRUE & (conf.int==TRUE | conf.band==TRUE)){
        lower.conf.int <- AUC-se*qnorm(1-(1-0.95)/2)
        upper.conf.int <- AUC+se*qnorm(1-(1-0.95)/2)
        lines(times,lower.conf.int,lty=2,lwd=2,col=col)
        lines(times,upper.conf.int,lty=2,lwd=2,col=col)
      }else{
        cat("Confidence intervals cannot be computed because you have chosen iid=FALSE (default) in the input of timeROC(). \n")
      }
    }
    if(conf.band==TRUE){
      if(object1$iid==TRUE & object2$iid==TRUE & (conf.int==TRUE | conf.band==TRUE)){
        vect.Delta <- rep(NA,2000)
        for (i in 1:2000){
          temp1 <- mat.iid*rnorm(object1$n)
          temp2 <- t(t(temp1)/se)
          vect.Delta[i] <- max(abs(colMeans(temp2)))   
        }
        C.alpha <- quantile(vect.Delta,0.95)
        lower.conf.band <- AUC-se*C.alpha
        upper.conf.band <- AUC+se*C.alpha
        lines(times,lower.conf.band,lty=3,lwd=2,col=col)
        lines(times,upper.conf.band,lty=3,lwd=2,col=col)
      }else{
        cat("Confidence band cannot be computed because you have chosen iid=FALSE (default) in the input of timeROC() for one of the two object. \n")
      }
    }
  }
  # }}}
  # {{{ class "ipcwcompetingrisksROC"
  if(class(object1)=="ipcwcompetingrisksROC"){
    if (FP==2){
      AUC <-  object1$AUC_2[!(is.na(object1$AUC_2)|is.na(object2$AUC_2))]- object2$AUC_2[!(is.na(object1$AUC_2)|is.na(object2$AUC_2))]
      times <- object1$times[!(is.na(object1$AUC_2)|is.na(object2$AUC_2))]
    }
    if (FP==1){
      AUC <-  object1$AUC_1[!(is.na(object1$AUC_1)|is.na(object2$AUC_1))]- object2$AUC_1[!(is.na(object1$AUC_1)|is.na(object2$AUC_1))]
      times <- object1$times[!(is.na(object1$AUC_1)|is.na(object2$AUC_1))]
    }
    if (FP!=1 & FP!=2){
      stop("FP indicates the type of specificity. It must be 1 or 2 \n")
    }
    if (add==FALSE){
      plot(times,AUC,xlab="time t",ylab=expression(paste(Delta,"AUC(t)")),lwd=2,ylim=ylim,type="l",col=col)
      abline(h=0,lty=2)
    }else{
      lines(times,AUC,lwd=2,type="l",col=col)
    }
    if(object1$iid==TRUE & object2$iid==TRUE & (conf.int==TRUE | conf.band==TRUE)){
      if (FP==1){
        mat.iid <- object1$inference$mat_iid_rep_1-object2$inference$mat_iid_rep_1
      }
      if (FP==2){
        mat.iid <- object1$inference$mat_iid_rep_2-object2$inference$mat_iid_rep_2
      }
      se<-apply(mat.iid,2,sd)/sqrt(object1$n)
    }
    if(conf.int==TRUE){
      if(object1$iid==TRUE & object2$iid==TRUE & (conf.int==TRUE | conf.band==TRUE)){
        lower.conf.int <- AUC-se*qnorm(1-(1-0.95)/2)
        upper.conf.int <- AUC+se*qnorm(1-(1-0.95)/2)
        lines(times,lower.conf.int,lty=2,lwd=2,col=col)
        lines(times,upper.conf.int,lty=2,lwd=2,col=col)
      }else{
        cat("Confidence band cannot be computed because you have chosen iid=FALSE (default) in the input of timeROC() for one of the two object. \n")
      }
    }
    if(conf.band==TRUE){
      if(object1$iid==TRUE & object2$iid==TRUE & (conf.int==TRUE | conf.band==TRUE)){
        vect.Delta <- rep(NA,2000)
        for (i in 1:2000){
          temp1 <- mat.iid*rnorm(object1$n)
          temp2 <- t(t(temp1)/se)
          vect.Delta[i] <- max(abs(colMeans(temp2)))   
        }
        C.alpha <- quantile(vect.Delta,0.95)
        lower.conf.band <- AUC-se*C.alpha
        upper.conf.band <- AUC+se*C.alpha
        lines(times,lower.conf.band,lty=3,lwd=2,col=col)
        lines(times,upper.conf.band,lty=3,lwd=2,col=col)
      }else{
        cat("Confidence band cannot be computed because you have chosen iid=FALSE (default) in the input of timeROC(). \n")
      }
    }
  }
  # }}}
  # {{{ error if not class ipcwcompetingrisksROC nor ipcwsurvivalROC
  if((class(object1)!="ipcwcompetingrisksROC" & class(object1)!="ipcwsurvivalROC") | (class(object2)!="ipcwcompetingrisksROC" & class(object2)!="ipcwsurvivalROC")){
    stop("Function written for object of class ipcwsurvivalROC or ipcwcompetingrisksROC\n")
  }
  # }}}
}
