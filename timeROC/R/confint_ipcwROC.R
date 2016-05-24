# {{{ input  :
# object : object od class confint.ipcwsurvivalROC or confint.ipcwcompetingrisksROC
# level : 1-alpha level
# n.sim : the number of Monte Carlo Simulation to estimate the quantile required to compute the confidence band
# }}}
# {{{ Output :
# Conf.band : matrix with rows containing upper and lower values of the simultaneous confidence band
# Conf.int : matrix with rows containing upper and lower values of the pointwise confidence interval
# }}}
confint.ipcwsurvivalROC <- function(object,parm=NULL,level=0.95,n.sim=2000,...){
  if(object$iid==FALSE){  
    stop(paste("Confidence intervals cannot be computed because you have chosen iid=FALSE (default) in the input of timeROC(). \n",sep=""))
  }
 ## browser()
  # {{{ remove NA if there is NA
  AUC <-  object$AUC[!is.na(object$AUC)]
  se <- object$inference$vect_sd_1[!is.na(object$AUC)]
  mat.iid <- object$inference$mat_iid_rep_1[,!is.na(object$AUC),drop=FALSE]
  # }}}
  # {{{ Pointwise confidence intervals
  lower <- AUC*100-se*100*qnorm(1-(1-level)/2)
  upper <- AUC*100+se*100*qnorm(1-(1-level)/2)
  tab_AUC_1<-round(cbind(lower,upper),2)
  colnames(tab_AUC_1)<-c(paste(((1-level)/2)*100,"%",sep=""),paste((1-(1-level)/2)*100,"%",sep=""))
  # }}}
  # {{{ Simultaneous confidence band
  # compute the quantile
  vect.Delta <- rep(NA,n.sim)
  for (i in 1:n.sim){
    temp1 <- mat.iid*rnorm(object$n)
    temp2 <- t(t(temp1)/se)
    vect.Delta[i] <- max(abs(colMeans(temp2)))   
  }
  C.alpha <- quantile(vect.Delta,level)
  # compute lower and upper bound
  lower.conf.band <- AUC*100-se*100*C.alpha
  upper.conf.band <- AUC*100+se*100*C.alpha
  ConfBand_AUC_1<-round(cbind(lower.conf.band,upper.conf.band),2)
  colnames(ConfBand_AUC_1)<-c(paste(((1-level)/2)*100,"%",sep=""),paste((1-(1-level)/2)*100,"%",sep=""))
  # }}}
  return(list(CI_AUC=tab_AUC_1,
              CB_AUC=ConfBand_AUC_1,
              C.alpha=C.alpha
              ))
}


confint.ipcwcompetingrisksROC <- function(object,parm=NULL,level=0.95,n.sim=2000,...){
  if(object$iid==FALSE){  
    stop(paste("Confidence intervals cannot be computed because you have chosen iid=FALSE (default) in the input of timeROC(). \n",sep=""))
  }
 ## browser()
  # {{{ remove NA if there is NA
  AUC.1 <-  object$AUC_1[!is.na(object$AUC_1)]
  se.1 <- object$inference$vect_sd_1[!is.na(object$AUC_1)]
  mat.iid.1 <- object$inference$mat_iid_rep_1[,!is.na(object$AUC_1),drop=FALSE]
  AUC.2 <-  object$AUC_2[!is.na(object$AUC_2)]
  se.2 <- object$inference$vect_sd_2[!is.na(object$AUC_2)]
  mat.iid.2 <- object$inference$mat_iid_rep_2[,!is.na(object$AUC_2),drop=FALSE]
  # }}}
  # {{{ Pointwise confidence intervals  
  lower_1<-AUC.1*100-se.1*100*qnorm(1-(1-level)/2)
  upper_1<-AUC.1*100+se.1*100*qnorm(1-(1-level)/2)
  tab_AUC_1<-round(cbind(lower_1,upper_1),2)
  colnames(tab_AUC_1)<-c(paste(((1-level)/2)*100,"%",sep=""),paste((1-(1-level)/2)*100,"%",sep=""))
  lower_2<-AUC.2*100-se.2*100*qnorm(1-(1-level)/2)
  upper_2<-AUC.2*100+se.2*100*qnorm(1-(1-level)/2)
  tab_AUC_2<-round(cbind(lower_2,upper_2),2)
  colnames(tab_AUC_2)<-c(paste(((1-level)/2)*100,"%",sep=""),paste((1-(1-level)/2)*100,"%",sep=""))
  # }}}
  # {{{ Simultaneous confidence band
  # compute the quantile
  vect.Delta.1 <- rep(NA,n.sim)
  vect.Delta.2 <- rep(NA,n.sim)
  for (i in 1:n.sim){
    rand.vect.norm <- rnorm(object$n)
    temp1.1 <- mat.iid.1*rand.vect.norm
    temp2.1 <- t(t(temp1.1)/se.1)
    vect.Delta.1[i] <- max(abs(colMeans(temp2.1)))
    temp1.2 <- mat.iid.2*rand.vect.norm
    temp2.2 <- t(t(temp1.2)/se.2)
    vect.Delta.2[i] <- max(abs(colMeans(temp2.2)))   
  }
  C.alpha.1 <- quantile(vect.Delta.1,level)
  C.alpha.2 <- quantile(vect.Delta.2,level)
  # compute lower and upper bound
  lower.conf.band.1 <- AUC.1*100-se.1*100*C.alpha.1
  upper.conf.band.1 <- AUC.1*100+se.1*100*C.alpha.1
  ConfBand_AUC_1<-round(cbind(lower.conf.band.1,upper.conf.band.1),2)
  colnames(ConfBand_AUC_1)<-c(paste(((1-level)/2)*100,"%",sep=""),paste((1-(1-level)/2)*100,"%",sep=""))
  lower.conf.band.2 <- AUC.2*100-se.2*100*C.alpha.2
  upper.conf.band.2 <- AUC.2*100+se.2*100*C.alpha.2
  ConfBand_AUC_2<-round(cbind(lower.conf.band.2,upper.conf.band.2),2)
  colnames(ConfBand_AUC_2)<-c(paste(((1-level)/2)*100,"%",sep=""),paste((1-(1-level)/2)*100,"%",sep=""))
  # }}}
  return(list(CI_AUC_1=tab_AUC_1,
              CB_AUC_1=ConfBand_AUC_1,
              C.alpha.1=C.alpha.1,
              CI_AUC_2=tab_AUC_2,
              CB_AUC_2=ConfBand_AUC_2,
              C.alpha.2=C.alpha.2)
         )
}
