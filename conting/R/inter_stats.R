inter_stats <-
function(object,cutoff=0.75,n.burnin=0,thin=1,prob.level=0.95){

if(cutoff<0 | cutoff>1){
stop("cutoff is a probability and should be between 0 and 1")}
if(n.burnin<0){
stop("n.burnin should be positive")}
if(n.burnin>=length(object$MODEL)){
stop("n.burnin should be less than the MCMC sample size (n.sample)")}
if(thin<1){
stop("thin should be greater than or equal to 1")}
if(prob.level<0 | prob.level>1){
stop("prob.level is a probability and should be between 0 and 1")}

if(n.burnin>0){
innerBETA<-object$BETA[-(1:n.burnin),]
innerMODEL<-object$MODEL[-(1:n.burnin)]} else{
innerBETA<-object$BETA
innerMODEL<-object$MODEL}

n.sample<-length(innerMODEL)
every<-seq(from=thin,to=n.sample,by=thin)
innerBETA<-innerBETA[every,]
innerMODEL<-innerMODEL[every]

INDO<-model2index(innerMODEL,dig=dim(innerBETA)[2])
innerBETA2<-innerBETA
innerBETA2[INDO==0]<-NA

post_mean<-apply(innerBETA2,2,mean,na.rm=TRUE)
post_var<-apply(innerBETA2,2,var,na.rm=TRUE)
post_prob<-apply(INDO,2,mean)
post_int<-c()
for(j in 1:length(post_mean)){
if(length(innerBETA[INDO[,j]==1,j])>1){
into<-HPDinterval(mcmc(innerBETA[INDO[,j]==1,j]),prob=prob.level)} else{
into<-c(NA,NA)}
post_int<-rbind(post_int,into)}

dimnames(post_int)<-NULL
lower<-post_int[,1]
upper<-post_int[,2]

est<-list(term=(dimnames(object$maximal.mod$x)[[2]])[post_prob>=cutoff],prob=post_prob[post_prob>=cutoff],
post_mean=post_mean[post_prob>=cutoff],post_var=post_var[post_prob>=cutoff],
lower=lower[post_prob>=cutoff],upper=upper[post_prob>=cutoff],thin=thin,prob.level=prob.level)

class(est)<-"interstat"

est}
