`coef.FitFGN` <-
function(object, ...){
BETA<-c(object$H, object$muHat)
sdB<-c(object$SEH,object$SEmu)
Z<-round((BETA-c(0.5,0))/sdB,3)
rn<-c("H", "mu")
cn<-c("MLE","sd","Z-ratio")
ans<-matrix(c(BETA,sdB, Z),ncol=3)
dimnames(ans)<-list(rn, cn)
ans
}

