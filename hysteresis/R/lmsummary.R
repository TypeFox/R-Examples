
lmsummary <-
function (g,N,studentize,center,cbb,joint) {
    wr1<- g$x-g$pred.x
    wr2<- g$y-g$pred.y
    if (center==TRUE) {
    wr1 <- wr1 - mean(wr1) 
    wr2 <- wr2 - mean(wr2)
    }
    n <- g$fit.statistics["n"]
    if (is.numeric(cbb)==TRUE){
      k <- n/cbb
      if (abs(k-round(k)) > 0.00001 || cbb <= 0 || abs(cbb-round(cbb)) > 0.00001) stop("invalid value for cbb.")}
    if (studentize==TRUE) {
      Xmat <- cbind(rep(1,n),sin(g$period.time),cos(g$period.time))
      h <- Xmat%*%solve(crossprod(Xmat))%*%t(Xmat)
      r.Ta <- wr1/sqrt(1-diag(h))
      r.Tb <- wr2/sqrt(1-diag(h))     
      wr1 <- r.Ta-mean(r.Ta)
      wr2 <- r.Tb-mean(r.Tb) 
    }
    hystenv$warningcountHysteresis <- 0

    bootdat <- mapply(lmbootwrapper, j=1:N, MoreArgs=list(wr1=wr1,wr2=wr2,x.pred=g$pred.x,y.pred=g$pred.y,n=n,cbb=cbb,joint=joint))
    bootint2 <- matrix(internal.1(bootdat[4,],bootdat[5,],bootdat[3,]),nrow=N,ncol=3)
    bootderived <- matrix(derived.1(bootdat[4,],bootdat[5,],bootdat[3,],bootint2[,1],bootint2[,2],bootint2[,3],rep(g$fit.statistics["period"],N)),nrow=N,ncol=3)
    bootamps <- matrix(derived.amps(bootint2[,1],bootint2[,2],bootint2[,3]),nrow=N,ncol=5)  
    bootfocus<- matrix(derived.focus(bootdat[4,],bootdat[5,],bootdat[3,]),nrow=N,ncol=3)
    bootdat <- data.frame(t(bootdat),bootderived,bootint2,bootamps,bootfocus,rep(g$fit.statistics["n"],N))
    colnames(bootdat) <- names(g$values)
    
if (diff(range(bootdat[,"rote.deg"])) > 170)
  warning("Bootstrapped rote.deg values on both sides of 0, 180 degrees.")
if (hystenv$warningcountHysteresis > 0) warning("Model failed to run ",hystenv$warningcountHysteresis," times.")
    
error<-apply(bootdat,2,sd,na.rm=T)
themean<-apply(bootdat,2,mean,na.rm=T)
ranges<-apply(bootdat,2,quantile,probs=c(0.025,0.25,0.5,0.75,0.975),na.rm=T)
full <- data.frame(g$values,t(ranges),error,themean)
colnames(full) <- c("Orig.Estimate","B.q0.025","B.q0.25","B.q0.5","B.q0.75","B.q0.975","Std.Error","Boot.Mean") 
    
full$Bias <- full$Boot.Mean-full$Orig.Estimate
full$Boot.Estimate <- full$Orig.Estimate-full$Bias  
    full[,c("B.q0.025","B.q0.25","B.q0.5","B.q0.75","B.q0.975")]<-full[,c("B.q0.025","B.q0.25","B.q0.5","B.q0.75","B.q0.975")]-
      2*full$Bias
    
    
    rad<-g$period.time
    pred.x<-full["b.x","Boot.Estimate"]*cos(rad)+full["cx","Boot.Estimate"] 
    pred.y<-full["b.y","Boot.Estimate"]*cos(rad)+full["retention","Boot.Estimate"]*sin(rad)+full["cy","Boot.Estimate"]
   
    full <- full[c("b.x","b.y",
                   "cx","cy","retention","coercion","area",
                   "lag","split.angle","hysteresis.x","hysteresis.y","ampx","ampy","rote.deg","rote.rad",
                   "semi.major","semi.minor","focus.x","focus.y","eccentricity"),]
    bootEst<-full[,"Boot.Estimate"]
    names(bootEst) <- rownames(full)
    bootStd<-full[,"Std.Error"]
    names(bootStd) <- rownames(full)
full2<-list("values"=full,"Boot.Estimates"=bootEst,"Boot.Std.Errors"=bootStd,
"fit.statistics"=g$fit.statistics,"pred.x"=g$pred.x,"pred.y"=g$pred.y, "x"=g$x,"y"=g$y,"call"=g$call,"method"=g$method,"boot.data"=bootdat,
            "Delta.Std.Errors"=g$Std.Errors)
invisible(full2)

}
