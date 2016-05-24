beta_mode <-
function(X,prior="SBH",y,IP,a=0.001,b=0.001){

Xt<-t(X)							## transpose of X
sy<-log(ifelse(y>0,y,1/6))					## starting values
sb<-coef(lm(sy~X-1))						## ""

###################### Log likelihoods ########################

loglik<-function(beta){
eta<-as.vector(X%*%matrix(beta,ncol=1))
sum(dpois(x=y,lambda=exp(eta),log=TRUE))}			## log-likelihood

dloglik<-function(beta){
eta<-as.vector(X%*%matrix(beta,ncol=1))
as.vector(Xt%*%matrix(y-exp(eta),ncol=1))}			## gradient of log-likelihood

d2loglik<-function(beta){
eta<-as.vector(X%*%matrix(beta,ncol=1))
w<-exp(eta)
-crossprod(x=X*w,y=X)}						## hessian of log-likelihood

################################################################

priortypes<-c("UIP","SBH")
priornum<-c(1,2)[prior==priortypes]				## which prior - defines prior part of posterior

################# UIP Priors ###################################

if(priornum==1){
iSig<-IP[-1,-1]				## inverse prior variance
#Sig<-solve(iSig)						## prior variance
Sig<-chol2inv(chol(iSig))						## prior variance
pp<-dim(X)[2]-1							## number of parameters - 1

prior<-function(beta){
dmvnorm(x=beta[-1],mean=rep(0,pp),sigma=Sig,log=TRUE)}		## prior log pdf

dprior<-function(beta){
-c(0,as.vector(iSig%*%matrix(beta[-1],ncol=1)))}		## gradient of prior log pdf

d2prior<-function(beta){
-cbind(0,rbind(0,iSig))}}					## hessian of prior log pdf

################# SBH Priors ###################################

if(priornum==2){						## parameters from inverse gamma
#a<-0.001
#b<-0.001
iSig<-IP[-1,-1]				## inverse scale matrix for t-distribution
pp<-dim(X)[2]-1

prior<-function(beta){
bsb<-as.vector(matrix(beta[-1],nrow=1)%*%iSig%*%matrix(beta[-1],ncol=1))
-0.5*(a+pp)*log(b+bsb)}						## prior log pdf				

dprior<-function(beta){
sb<-as.vector(iSig%*%matrix(beta[-1],ncol=1))
bsb<-sum(beta[-1]*sb)
-c(0,(a+pp)*sb/(b+bsb))}					## gradient of prior log pdf

d2prior<-function(beta){
sb<-as.vector(iSig%*%matrix(beta[-1],ncol=1))
bsb<-sum(beta[-1]*sb)						## hessian of prior log pdf
-cbind(0,rbind(0,(a+pp)*((iSig/(b+bsb))-2*((matrix(sb,ncol=1)%*%matrix(sb,nrow=1))/((b+bsb)^2)))))}}

##################################################################

mlogpost<-function(beta){
-loglik(beta)-prior(beta)}					## minus log posterior

mdlogpost<-function(beta){
-dloglik(beta)-dprior(beta)}					## gradient of minus log posterior

md2logpost<-function(beta){
-d2loglik(beta)-d2prior(beta)}					## hessian of minus log posterior

opt<-nlminb(start=sb,objective=mlogpost,gradient=mdlogpost,hessian=md2logpost) ## optimisation!

mbeta<-opt$par							

mbeta}
