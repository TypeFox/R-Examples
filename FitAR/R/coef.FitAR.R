"coef.FitAR" <-
function(object, ...){
if (object$SubsetQ)
    if (object$FitMethod=="MLE")
        BETA<-object$zetaHat
    else
        BETA<-(object$phiHat)[object$pvec]
else
    BETA<-object$phiHat
p<-length(BETA)
BETA<-c(BETA,object$muHat)
sdB<-sqrt(diag(object$covHat))
sdB<-c(sdB,sqrt((object$sigsq)/(length(object$res)*sum(c(1,-object$phiHat))^2)))
Z<-BETA/sdB
if (object$SubsetQ)
    if (object$FitMethod=="MLE")
        rn<-c(paste("zeta(",object$pvec,")", sep=""),"mu")
    else
        rn<-c(paste("phi(",object$pvec,")", sep=""),"mu")
    
else
    rn<-c(paste("phi(",1:p,")", sep=""),"mu")
cn<-c("MLE","sd","Z-ratio")
ans<-matrix(c(BETA,sdB, Z),ncol=3)
dimnames(ans)<-list(rn, cn)
ans
}
