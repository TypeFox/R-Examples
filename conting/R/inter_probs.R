inter_probs <-
function(object,cutoff=0.75,n.burnin=0,thin=1){

if(n.burnin<0){
stop("n.burnin should be positive")}
if(n.burnin>=length(object$MODEL)){
stop("n.burnin should be less than the MCMC sample size (n.sample)")}
if(thin<1){
stop("thin should be greater than or equal to 1")}
if(cutoff<0 | cutoff>1){
stop("cutoff is a probability and should be between 0 and 1")}

term.labels<-c("(Intercept)",attr(summary(object$maximal.mod)$terms,"term.labels"))
term.numbers<-attributes(object$maximal.mod$x)$assign
small<-function(i){
min((1:length(term.numbers))[term.numbers==i])}
ids<-sapply(X=0:max(term.numbers),FUN=small)

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

post_prob<-apply(INDO[,ids],2,mean)

est<-list(term=term.labels[post_prob>=cutoff],prob=post_prob[post_prob>=cutoff],thin=thin)

class(est)<-"interprob"

est

}
