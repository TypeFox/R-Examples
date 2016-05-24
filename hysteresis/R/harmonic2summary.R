harmonic2summary <-
function (g,N,studentize,cbb,joint) {

xresid<-g$pred.x-g$x
yresid<-g$pred.y-g$y
n <- length(xresid)-3
if (is.numeric(cbb)==TRUE){
k <- n/cbb
if (abs(k-round(k)) > 0.00001) stop("number of observations - 3 divided by cbb needs to be an integer.")}
if (studentize==TRUE) {
  h.Ta <- lm.influence(g$fit[[1]])$hat
  h.Tb <- lm.influence(g$fit[[2]])$hat
  r.Ta <- xresid/sqrt(1-h.Ta)
  r.Tb <- yresid/sqrt(1-h.Tb)     
  xresid <- r.Ta-mean(r.Ta)
  yresid <- r.Tb-mean(r.Tb) 
}

bootdat<-mapply(harmonic2boot,j=1:N,MoreArgs=list(pred.x=g$pred.x,pred.y=g$pred.y,xresid=xresid,yresid=yresid,ti=g$period.time,n=n,cbb=cbb,joint=joint))
bootdat<-t(bootdat)
bootderived <- matrix(derived.2(bootdat[,3],bootdat[,4],bootdat[,6],rep(g$fit.statistics["period"],N)),nrow=N,ncol=3)
bootinternal1 <- matrix(internal.2(bootdat[,3],bootdat[,4],bootdat[,6],bootdat[,5]),nrow=N,ncol=4)
bootamps <- matrix(derived.amps(bootdat[,3],bootdat[,4],bootdat[,6]),nrow=N,ncol=5)  
bootfocus <- matrix(derived.focus(bootinternal1[,3],bootinternal1[,4],bootinternal1[,1]),nrow=N,ncol=3)
bootdat <- data.frame(bootdat,bootderived,bootinternal1,bootamps,bootfocus)  
colnames(bootdat) <- names(g$values)

themean<-apply(bootdat,2,mean,na.rm=TRUE)

error<-apply(bootdat,2,sd,na.rm=TRUE)
ranges<-apply(bootdat,2,quantile,probs=c(0.025,0.25,0.5,0.75,0.975),na.rm=TRUE)
full <- data.frame(g$values,t(ranges),error, themean)
colnames(full) <- c("Orig.Estimate","B.q0.025","B.q0.25","B.q0.5","B.q0.75","B.q0.975","Std.Error","Boot.Mean")

full$Bias <- full$Boot.Mean-full$Orig.Estimate
full$Boot.Estimate <- full$Orig.Estimate-full$Bias
full[,c("B.q0.025","B.q0.25","B.q0.5","B.q0.75","B.q0.975")]<-full[,c("B.q0.025","B.q0.25","B.q0.5","B.q0.75","B.q0.975")]-
  2*full$Bias


if (diff(range(bootdat[,"rote.deg"])) > 170)
  warning("Bootstrapped rote.deg values on both sides of 0, 180 degrees.")
  

rad<-g$period.time+full["phase.angle","Boot.Estimate"]
pred.x<-full["b.x","Boot.Estimate"]*cos(rad)+full["cx","Boot.Estimate"] 
pred.y<-full["b.y","Boot.Estimate"]*cos(rad)+full["retention","Boot.Estimate"]*sin(rad)+full["cy","Boot.Estimate"]


full <- full[c("b.x","b.y","phase.angle",
               "cx","cy","retention","coercion","area",
               "lag","split.angle","hysteresis.x","hysteresis.y","ampx","ampy","rote.deg","rote.rad",
               "semi.major","semi.minor","focus.x","focus.y","eccentricity"),]
bootEst<-full[,"Boot.Estimate"]
names(bootEst) <- rownames(full)
bootStd<-full[,"Std.Error"]
names(bootStd) <- rownames(full)
full2<-list("values"=full,"Boot.Estimates"=bootEst,"Boot.Std.Errors"=bootStd,
"fit.statistics"=g$fit.statistics, "pred.x"=pred.x,"pred.y"=pred.y, "x"=g$x,"y"=g$y,"call"=g$call,"method"=g$method,"boot.data"=bootdat,
            "Delta.Std.Errors"=g$Std.Errors)
invisible(full2)
}
