total_pop <-
function(object,n.burnin=0,thin=1,prob.level=0.95){

if(n.burnin<0){
stop("n.burnin should be positive")}
if(n.burnin>=length(object$MODEL)){
stop("n.burnin should be less than the MCMC sample size (n.sample)")}
if(thin<1){
stop("thin should be greater than or equal to 1")}
if(prob.level<0 | prob.level>1){
stop("prob.level is a probability and should be between 0 and 1")}

missing<-c(object$missing1,object$missing2)

obs.z<-sum(object$maximal.mod$y[-missing])

if(is.matrix(object$Y0)){
if(n.burnin>0){
YY0<-matrix(object$Y0[-(1:n.burnin),],ncol=dim(object$Y0)[2])} else{
YY0<-object$Y0}

TOT<-apply(YY0,1,sum)+obs.z} else{

if(n.burnin>0){
YY0<-object$Y0[-(1:n.burnin)]} else{
YY0<-object$Y0}

TOT<-YY0+obs.z}

n.sample<-length(TOT)
every<-seq(from=thin,to=n.sample,by=thin)
TOT<-TOT[every]

int<-HPDinterval(mcmc(TOT),prob=prob.level)

est<-list(TOT=TOT,int=int,meanTOT=mean(TOT),thin=thin,prob.level=prob.level)

class(est)<-"totpop"

est}
