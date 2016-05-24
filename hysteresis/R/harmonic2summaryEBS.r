#Experimental
# harmonic2summaryEBS <-
# function (g,N,studentize,cbb,joint,EBS) {
# 
# xresid<-g$pred.x-g$x
# yresid<-g$pred.y-g$y
# n <- length(xresid)-3
# if (is.numeric(cbb)==TRUE){
# k <- n/cbb
# if (abs(k-round(k)) > 0.00001) stop("number of observations - 3 divided by cbb needs to be an integer.")}
# if (studentize==TRUE) {
#   h.Ta <- lm.influence(g$fit[[1]])$hat
#   h.Tb <- lm.influence(g$fit[[2]])$hat
#   r.Ta <- xresid/sqrt(1-h.Ta)
#   r.Tb <- yresid/sqrt(1-h.Tb)     
#   xresid <- r.Ta-mean(r.Ta)
#   yresid <- r.Tb-mean(r.Tb) 
# }
# 
# bootdat<-mapply(harmonic2boot,j=1:N,MoreArgs=list(pred.x=g$pred.x,pred.y=g$pred.y,xresid=xresid,yresid=yresid,ti=g$period.time,n=n,cbb=cbb,joint=joint))
# bootdat<-t(bootdat)
# bootderived <- matrix(derived.2(bootdat[,3],bootdat[,4],bootdat[,6],rep(g$fit.statistics["period"],N)),nrow=N,ncol=3)
# bootinternal1 <- matrix(internal.2(bootdat[,3],bootdat[,4],bootdat[,6],bootdat[,5]),nrow=N,ncol=4)
# bootamps <- matrix(derived.amps(bootdat[,3],bootdat[,4],bootdat[,6]),nrow=N,ncol=5)  
# bootfocus <- matrix(derived.focus(bootinternal1[,3],bootinternal1[,4],bootinternal1[,1]),nrow=N,ncol=3)
# bootdat <- data.frame(bootdat,bootderived,bootinternal1,bootamps,bootfocus)  
# colnames(bootdat) <- names(g$values)
# 
# themean<-apply(bootdat,2,mean,na.rm=TRUE)
# 
# error<-apply(bootdat,2,sd,na.rm=TRUE)
# ranges<-apply(bootdat,2,quantile,probs=c(0.025,0.25,0.5,0.75,0.975),na.rm=TRUE)
# full <- data.frame(g$values,t(ranges),error, themean)
# colnames(full) <- c("Orig.Estimate","B.q0.025","B.q0.25","B.q0.5","B.q0.75","B.q0.975","Std.Error","Boot.Mean")
# 
# full$Bias <- full$Boot.Mean-full$Orig.Estimate
# full$Boot.Estimate <- full$Orig.Estimate-full$Bias
# full[,c("B.q0.025","B.q0.25","B.q0.5","B.q0.75","B.q0.975")]<-full[,c("B.q0.025","B.q0.25","B.q0.5","B.q0.75","B.q0.975")]-
#   2*full$Bias
# 
# if (EBS==TRUE) {
# alpha <- full$Boot.Estimate-3*full$Std.Error
# beta <- full$Boot.Estimate+3*full$Std.Error
# ranparams<- matrix(runif(length(alpha)*200,alpha,beta),length(alpha),200)
# ranparams[7:9,] <- t(matrix(derived.2(ranparams[3,],ranparams[4,],ranparams[6,],rep(g$fit.statistics["period"],200)),nrow=200,ncol=3))
# ranparams[10:13,] <- t(matrix(internal.2(ranparams[3,],ranparams[4,],ranparams[6,],ranparams[5,]),nrow=200,ncol=4))
# ranparams[14:18,] <- t(matrix(derived.amps(ranparams[3,],ranparams[4,],ranparams[6,]),nrow=200,ncol=5))  
# ranparams[19:21,] <- t(matrix(derived.focus(ranparams[12,],ranparams[13,],ranparams[10,]),nrow=200,ncol=3))
# biasmat<- matrix(NA,length(alpha),200)
# themeanmat<- matrix(NA,length(alpha),200)
# for (i in 1:200) {
# themean2 <- ranparams[,i]
# rad2<-g$period.time+themean2[5]
#       pred.x2 <- themean2[3]*cos(rad2)+themean2[1]
#       pred.y2 <- themean2[4]*cos(rad2)+themean2[6]*sin(rad2)+themean2[2]
#      xresid2<-pred.x2-g$x
#       yresid2<-pred.y2-g$y
#       n <- length(xresid)-3
#   if (is.numeric(cbb)==TRUE){
# k <- n/cbb
# if (abs(k-round(k)) > 0.00001) stop("number of observations - 3 divided by cbb needs to be an integer.")}
# if (studentize==TRUE) {
#   h.Ta <- lm.influence(g$fit[[1]])$hat
#   h.Tb <- lm.influence(g$fit[[2]])$hat
#   r.Ta <- xresid/sqrt(1-h.Ta)
#   r.Tb <- yresid/sqrt(1-h.Tb)     
#   xresid2 <- r.Ta-mean(r.Ta)
#   yresid2 <- r.Tb-mean(r.Tb) 
# }
#     
#  bootdatrep<-mapply(harmonic2boot,j=1:100,MoreArgs=list(pred.x=pred.x2,pred.y=pred.y2,xresid=xresid2,yresid=yresid2,ti=g$period.time,n=n,cbb=cbb,joint=joint))
#  bootdatrep<-t(bootdatrep)
#  bootderivedrep <- matrix(derived.2(bootdatrep[,3],bootdatrep[,4],bootdatrep[,6],rep(g$fit.statistics["period"],100)),nrow=100,ncol=3)
#  bootinternal1rep <- matrix(internal.2(bootdatrep[,3],bootdatrep[,4],bootdatrep[,6],bootdatrep[,5]),nrow=100,ncol=4)
#  bootampsrep <- matrix(derived.amps(bootdatrep[,3],bootdatrep[,4],bootdatrep[,6]),nrow=100,ncol=5)  
#  bootfocusrep <- matrix(derived.focus(bootinternal1rep[,3],bootinternal1rep[,4],bootinternal1rep[,1]),nrow=100,ncol=3)
#  bootdatrep <- data.frame(bootdatrep,bootderivedrep,bootinternal1rep,bootampsrep,bootfocusrep)  
#  colnames(bootdatrep) <- names(g$values)
#  themeanrep<-apply(bootdatrep,2,mean,na.rm=TRUE)
#  themeanmat[,i] <- themeanrep
#  biasmat[,i] <- themeanrep-themean2
# }
# bmcoefs <- matrix(NA,nrow=dim(biasmat)[1],ncol=dim(biasmat)[1]+1)
# for (j in 1:dim(biasmat)[1]) {
# biasmodel <- try(lars(x=cbind(rep(1,200),t(themeanmat)),y=as.vector(biasmat[j,]))) 
# try(bmcoefs[j,] <- coef(biasmodel)[which.min(summary(biasmodel)$Cp),])
# }
# bmcoefs[is.na(bmcoefs)] <- 0
# OREst2 <- full$Orig.Estimate
# OREst2[is.na(OREst2)] <- 0
# full$Bias <- bmcoefs%*%c(1,full$Orig.Estimate)
# full$Boot.Estimate <- full$Orig.Estimate-full$Bias
# full[,c("B.q0.025","B.q0.25","B.q0.5","B.q0.75","B.q0.975")]<-full[,c("B.q0.025","B.q0.25","B.q0.5","B.q0.75","B.q0.975")]-
#   2*full$Bias
# }
# 
# if (diff(range(bootdat[,"rote.deg"])) > 170)
#   warning("Bootstrapped rote.deg values on both sides of 0, 180 degrees.")
#   
# 
# rad<-g$period.time+full["phase.angle","Boot.Estimate"]
# pred.x<-full["b.x","Boot.Estimate"]*cos(rad)+full["cx","Boot.Estimate"] 
# pred.y<-full["b.y","Boot.Estimate"]*cos(rad)+full["retention","Boot.Estimate"]*sin(rad)+full["cy","Boot.Estimate"]
# 
# 
# full <- full[c("b.x","b.y","phase.angle",
#                "cx","cy","retention","coercion","area",
#                "lag","split.angle","hysteresis.x","hysteresis.y","ampx","ampy","rote.deg","rote.rad",
#                "semi.major","semi.minor","focus.x","focus.y","eccentricity"),]
# bootEst<-full[,"Boot.Estimate"]
# names(bootEst) <- rownames(full)
# bootStd<-full[,"Std.Error"]
# names(bootStd) <- rownames(full)
# full2<-list("values"=full,"Boot.Estimates"=bootEst,"Boot.Std.Errors"=bootStd,
# "fit.statistics"=g$fit.statistics, "pred.x"=pred.x,"pred.y"=pred.y, "x"=g$x,"y"=g$y,"call"=g$call,"method"=g$method,"boot.data"=bootdat,
#             "Delta.Std.Errors"=g$Std.Errors)
# invisible(full2)
# }
