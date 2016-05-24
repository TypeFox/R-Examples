summary.bcct <-
function(object,n.burnin=0,thin=1,cutoff=0.75,statistic="X2",best=NULL,scale=0.1,prob.level=0.95,...){

if(n.burnin<0){
stop("n.burnin should be positive")}
if(n.burnin>=length(object$MODEL)){
stop("n.burnin should be less than the MCMC sample size (n.sample)")}
if(thin<1){
stop("thin should be greater than or equal to 1")}
if(cutoff<0 | cutoff>1){
stop("cutoff is a probability and should be between 0 and 1")}
if(statistic!="X2" & statistic!="FreemanTukey" & statistic!="deviance"){
stop("statistic not found")}
if(scale<0 | scale>1){
stop("scale should be between 0 and 1")}
if(!is.null(best)){
if(best<=0){
stop("best should be positive")}}
if(prob.level<0 | prob.level>1){
stop("prob.level is a probability and should be between 0 and 1")}

is1<-inter_stats(object,n.burnin=n.burnin,cutoff=cutoff,thin=thin,prob.level=prob.level)
is2<-mod_probs(object,n.burnin=n.burnin,scale=scale,best=best,thin=thin)
is3<-bayespval(object,n.burnin=n.burnin,thin=thin,statistic=statistic)

est<-list(BETA=object$BETA,MODEL=object$MODEL,SIG=object$SIG,rj_acc=object$rj_acc,mh_acc=object$mh_acc,priornum=object$priornum,maximal.mod=object$maximal.mod,IP=object$IP,eta.hat=object$eta.hat,save=object$save,name=object$name,int_stats=is1,mod_stats=is2,pval_stats=is3)

class(est)<-"sbcct"

est

}
