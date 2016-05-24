summary.loop2r <- function(object,boot=TRUE,N=1000,cbb=NULL,joint=FALSE,seed=NULL,...) {
  g <- object
  summarycall <- match.call()
  if (!is.null(seed)) set.seed(seed)
  if (boot==TRUE) {
    xresid<-g$pred.x-g$x
    yresid<-g$pred.y-g$y
    obs <- length(xresid)-3
    if (is.numeric(cbb)==TRUE){
      k <- obs/cbb
      if (abs(k-round(k)) > 0.00001) stop("number of observations - 3 divided by cbb needs to be an integer.")}
    
    if (g$method=="harmonic2")
    bootdat<-mapply(splitloopboot,j=1:N,MoreArgs=list(pred.x=g$pred.x,pred.y=g$pred.y,xresid=xresid,yresid=yresid,ti=g$period.time,obs=obs,n=g$values["n"],m=g$values["m"],extended.classical=g$extended.classical,cbb=cbb,joint=joint,period=g$period))
    else bootdat<-mapply(splitloopbootgeom,j=1:N,MoreArgs=list(pred.x=g$pred.x,pred.y=g$pred.y,xresid=xresid,yresid=yresid,obs=obs,extended.classical=g$extended.classical,cbb=cbb,joint=joint,period=g$period))

    bootdat<-t(bootdat)
    colnames(bootdat) <- names(g$values)
    
    error<-apply(bootdat,2,sd,na.rm=T)
    error <- ifelse(is.na(g$values),NA,error)
    ranges<-apply(bootdat,2,quantile,probs=c(0.025,0.25,0.5,0.75,0.975),na.rm=T)
    themean<-apply(bootdat,2,mean,na.rm=T)
    full <- data.frame(g$values,t(ranges),error, themean,g$Std.Errors)
    colnames(full) <- c("Orig.Estimate","B.q0.025","B.q0.25","B.q0.5","B.q0.75","B.q0.975","Std.Error","Boot.Mean","Delta.Error")
    
    full$Bias <- full$Boot.Mean-full$Orig.Estimate
    full$Boot.Estimate <- full$Orig.Estimate-full$Bias
    full[,c("B.q0.025","B.q0.25","B.q0.5","B.q0.75","B.q0.975")]<-full[,c("B.q0.025","B.q0.25","B.q0.5","B.q0.75","B.q0.975")]-
      2*full$Bias
    
    rad<-g$period.time+full["phase.angle","Boot.Estimate"]
    Ind <- (rad > 0) & (rad < 0)
    pred.x<-full["b.x","Boot.Estimate"]*cos(rad)+full["cx","Boot.Estimate"] 
    if (g$extended.classical==FALSE & g$method=="harmonic2")
      pred.y<-full["b.y","Boot.Estimate"]*cos(rad)^g$values["n"]+Ind*full["retention.above","Boot.Estimate"]*sin(rad)^g$values["m"]+(1-Ind)*full["retention.below","Boot.Estimate"]*sin(rad)^g$values["m"]+full["cy","Boot.Estimate"]
    else
      pred.y<-sign(cos(rad))*full["b.y","Boot.Estimate"]*abs(cos(rad))^full["n","Boot.Estimate"]+sign(sin(rad))*Ind*full["retention.above","Boot.Estimate"]*abs(sin(rad))^full["m","Boot.Estimate"]+sign(sin(rad))*(1-Ind)*full["retention.below","Boot.Estimate"]*abs(sin(rad))^g$values["m"]+full["cy","Boot.Estimate"]
   bootEst<-full[,"Boot.Estimate"]
    names(bootEst) <- rownames(full)
    bootStd<-full[,"Std.Error"]
    names(bootStd) <- rownames(full)
    full2<-list("values"=full,"Boot.Estimates"=bootEst,"Boot.Std.Errors"=bootStd,
                "pred.x"=pred.x,"pred.y"=pred.y, "x"=g$x,"y"=g$y,"call"=g$call,"extended.classical"=g$extended.classical,"boot.data"=bootdat,"boot"=TRUE,"method"=g$method)
    full2$summarycall <- summarycall
    full2$residuals <- sqrt((g$x-pred.x)^2+(g$y-pred.y)^2)
    full2$period.time <- g$period.time
  class(full2) <- "loop2rsummary"  
    full2 
  }
  else { res <- g
         res$summarycall <- summarycall
  res}
}
