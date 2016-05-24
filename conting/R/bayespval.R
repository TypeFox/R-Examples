bayespval <-
function(object,n.burnin=0,thin=1,statistic="X2"){

if(n.burnin<0){
stop("n.burnin should be positive")}
if(n.burnin>=length(object$MODEL)){
stop("n.burnin should be less than the MCMC sample size (n.sample)")}
if(thin<1){
stop("thin should be greater than or equal to 1")}
if(statistic!="X2" & statistic!="FreemanTukey" & statistic!="deviance"){
stop("statistic not found")}

if(is.null(object$missing1)){
yyy<-object$maximal.mod$y
xxx<-object$maximal.mod$x} else{
missing<-c(object$missing1,object$missing2)
yyy<-object$maximal.mod$y[-missing]
xxx<-object$maximal.mod$x[-missing,]}

if(n.burnin>0){
innerBETA<-object$BETA[-(1:n.burnin),]} else{
innerBETA<-object$BETA}

n.sample<-dim(innerBETA)[1]

every<-seq(from=thin,to=n.sample,by=thin)

innerBETA<-innerBETA[every,]

MU<-exp(innerBETA%*%t(xxx))
PRED<-matrix(rpois(n=prod(dim(MU)),lambda=as.vector(MU)),ncol=dim(MU)[2])
Y<-matrix(rep(yyy,dim(MU)[1]),ncol=dim(MU)[2],byrow=TRUE)

statnum<-(1:3)[c("X2","deviance","FreemanTukey")==statistic]

if(statnum==1){
Tpred<-apply(((PRED-MU)^2)/MU,1,sum)
Tobs<-apply(((Y-MU)^2)/MU,1,sum)}

if(statnum==2){
Tpred<-apply(-2*dpois(x=PRED,lambda=MU,log=TRUE),1,sum)
Tobs<-apply(-2*dpois(x=Y,lambda=MU,log=TRUE),1,sum)}

if(statnum==3){
Tpred<-apply((sqrt(PRED)-sqrt(MU))^2,1,sum)
Tobs<-apply((sqrt(Y)-sqrt(MU))^2,1,sum)}

pval<-mean(as.numeric(Tpred>Tobs))

est<-list(PRED=PRED,Tpred=Tpred,Tobs=Tobs,pval=pval,statnum=statnum,statistic=statistic,thin=thin)

class(est)<-"pval"

est}
