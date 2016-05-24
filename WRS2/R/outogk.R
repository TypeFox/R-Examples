outogk<-function(x,sigmamu=taulc,v=gkcov,op=TRUE,SEED=FALSE,
beta=max(c(.95,min(c(.99,1/nrow(x)+.94)))),n.iter=1,plotit=TRUE,...){
#
# Use the ogk estimator to
# determine which points are outliers
#
#  op=T uses robust Mahalanobis distance based on
#  the OGK estimator with  beta adjusted so that
#  the outside rate per observation is approximately .05
#  under normality.
#  op=F returns the outliers based on the distances used
#  by the OGK estimator
#  (Currently, op=T seems best for detecting outliers.)
#
if(!is.matrix(x))stop("x should be a matrix")
x<-elimna(x)
if(!op){
temp<-ogk.pairwise(x,sigmamu=sigmamu,v=v,beta=beta,n.iter=n.iter,...)
vals<-hard.rejection(temp$distances,p=ncol(x),beta=beta,...)
flag<-(vals==1)
vals<-c(1:nrow(x))
outid<-vals[!flag]
keep<-vals[flag]
if(is.matrix(x)){
if(ncol(x)==2 && plotit){
plot(x[,1],x[,2],xlab="X", ylab="Y",type="n")
points(x[flag,1],x[flag,2])
if(sum(!flag)>0)points(x[!flag,1],x[!flag,2],pch="o")
}}}
if(op){
temp<-out(x,cov.fun=ogk,beta=beta,plotit=plotit,SEED=SEED)
outid<-temp$out.id
keep<-temp$keep
}
list(out.id=outid,keep=keep,distances=temp$dis)
}