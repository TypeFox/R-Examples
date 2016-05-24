

 ###################
 #### plot.dblm ####
 ###################

 ## Description:
 ##    generic function plot for a dblm object. Six plots (selected by which)
 ##    are available:
 ##       - a plot of residual vs fitted values.
 ##       - the Q-Qplot of normality.
 ##       - a Scale-Location plot of residuals against fitted values.
 ##       - cook.distance
 ##       - Leverage vs residuals
 ##        - minimum effective rank of OCV, GCV, AIC or BIC method.
 ##


plot.dbplsr<-function(x,which=c(1L:4L),main="",scores.comps=1:2,component=1,               
              method=c("OCV","GCV","AIC","BIC"),...){

     # stop if the object is not a dblm object.
     if (!inherits(x, "dbplsr"))
        stop("\"dbplsr\" objects")
        
      
     # auxiliar boolean vector variable with 4 items. True if the plot is
     # selected by which or false if not.
     if (!is.numeric(which) || any(which < 1) || any(which > 4))
        stop("'which' must be in 1:4")

     show <- rep(FALSE, 4)
     show[which] <- TRUE

     # change the panel like lm plot (with mouse clik or with enter).
     # if only one plot is selected ask=FALSE --> not need to change the panel.
     ask = prod(par("mfcol")) < length(which)
      one.fig <- prod(par("mfcol")) == 1
     if (ask) {
        oask <- devAskNewPage(TRUE)  # ask for new page
        on.exit(devAskNewPage(oask)) # exit the format plot (if true)
     }
       
     # pairs between components scores   
     if(show[1L]){
       if ((main==""&&length(which)==1)||length(which)>1)
        main<-"Score plot"
       plot(as.data.frame(x$fk)[,scores.comps],main=main)
     }
     
     # pairs between y and components scores  
     if(show[2L]){
       if ((main==""&&length(which)==1)||length(which)>1)
        main<-paste("Response vs score",component)
       xlab=paste("score ",component)
       ylab="response"     
       y<-x$y
       fk <- x$fk[[component]]
       plot(fk,y,main=main,xlab=xlab,ylab=ylab)
       dblmplot<-dblm(y~fk,method="eff.rank",eff.rank=1)
       aa<- sort(fk, decreasing = FALSE,index.return = T)
       fittord<-dblmplot$fitt[aa$ix]
       lines(aa$x,fittord,col="red")
     }
      
     # contribution to final R2     
     if(show[3L]){
      R2<-summary(x)$r.squared
      R2comp<-c(R2[1],diff(R2))
      if ((main==""&&length(which)==1)||length(which)>1)
        main<-"R2 contribution"
      plot(R2comp,xlab="Components",ylab="R squared", main=main,type="o")  
     } 
     
     # components vs method statistical 
     if(show[4L]){
      xlab="Components"
      if ((main==""&&length(which)==1)||length(which)>1)
        main<-paste(method[1], "vs Components ")
      if (method[1]=="OCV")
       plot(x$ocv,type="o",xlab=xlab,ylab="OCV",main=main)
      if (method[1]=="GCV")
       plot(x$gcv,type="o",xlab=xlab,ylab="GCV",main=main)
      if (method[1]=="AIC")
       plot(x$aic,type="o",xlab=xlab,ylab="AIC",main=main)
      if (method[1]=="BIC")
       plot(x$bic,type="o",xlab=xlab,ylab="BIC",main=main)
    }   
  return (invisible())
}