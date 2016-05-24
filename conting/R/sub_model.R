sub_model <-
function(object,formula=NULL,order=1,n.burnin=0,thin=1,prob.level=0.95,statistic="X2"){

if(n.burnin<0){
stop("n.burnin should be positive")}
if(n.burnin>=length(object$MODEL)){
stop("n.burnin should be less than the MCMC sample size (n.sample)")}
if(thin<1){
stop("thin should be greater than or equal to 1")}
if(prob.level<0 | prob.level>1){
stop("prob.level is a probability and should be between 0 and 1")}
if(order<1){
stop("order should be greater than or equal to 1")}
if(statistic!="X2" & statistic!="FreemanTukey" & statistic!="deviance"){
stop("statistic not found")}

if(n.burnin>0){
if(class(object)=="bict"){
innerY0<-matrix(object$Y0[-(1:n.burnin),],ncol=dim(object$Y0)[2])} else{
innerY0<-NULL}
innerSIG<-object$SIG[-(1:n.burnin)]
innerBETA<-object$BETA[-(1:n.burnin),]
innerMODEL<-object$MODEL[-(1:n.burnin)]} else{
if(class(object)=="bict"){
innerY0<-matrix(object$Y0,ncol=dim(object$Y0)[2])} else{
innerY0<-NULL}
innerSIG<-object$SIG
innerBETA<-object$BETA
innerMODEL<-object$MODEL}

n.sample<-length(innerMODEL)
every<-seq(from=thin,to=n.sample,by=thin)
innerBETA<-innerBETA[every,]
innerMODEL<-innerMODEL[every]
innerSIG<-innerSIG[every]
if(!is.null(innerY0)){
innerY0<-matrix(innerY0[every,],ncol=dim(object$Y0)[2])}

tab<-sort(table(innerMODEL)/length(innerMODEL),decreasing=TRUE)

if(is.null(formula)){
border<-order
if(order>length(tab)){
stop("Model not visited in (thinned) sample")}
interest<-tab[order]
bformula<-index2formula(index=model2index(model=names(interest),dig=dim(object$maximal.mod$x)[2]),maximal.mod=object$maximal.mod)} else{
bformula<-formula
interest<-tab[names(tab)==index2model(formula2index(big.X=object$maximal.mod$x, formula=formula, data=object$maximal.mod$data))]
if(length(interest)==0){
stop("Model not visited in (thinned) sample")}
border<-(1:length(tab))[names(tab)==index2model(formula2index(big.X=object$maximal.mod$x, formula=formula, data=object$maximal.mod$data))]
}
bformula<-paste0("~",as.character(bformula)[3])

int_index<-model2index(model=names(interest),dig=dim(object$maximal.mod$x)[2])

redBETA<-innerBETA[innerMODEL==names(interest),int_index==1]
redSIG<-innerSIG[innerMODEL==names(interest)]
if(!is.null(innerY0)){
redY0<-matrix(innerY0[innerMODEL==names(interest),],ncol=dim(object$Y0)[2])} else{
redY0<-NULL}

if(!is.null(redY0)){
missing<-c(object$missing1,object$missing2)
yyy<-object$maximal.mod$y[-missing]
xxx<-object$maximal.mod$x[-missing,int_index==1]
obs.z<-sum(yyy)
redTOT<-apply(redY0,2,sum)+obs.z
int<-HPDinterval(mcmc(redTOT),prob=prob.level)
meanTOT<-mean(redTOT)} else{
yyy<-object$maximal.mod$y
xxx<-object$maximal.mod$x[,int_index==1]
redTOT<-NULL
int<-NULL
meanTOT<-NULL}

MU<-exp(redBETA%*%t(xxx))
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

post_mean<-apply(redBETA,2,mean,na.rm=TRUE)
post_var<-apply(redBETA,2,var,na.rm=TRUE)
post_int<-c()
for(j in 1:length(post_mean)){
into<-HPDinterval(mcmc(redBETA[,j]),prob=prob.level)
post_int<-rbind(post_int,into)}

dimnames(post_int)<-NULL
lower<-post_int[,1]
upper<-post_int[,2]

est<-list(term=(dimnames(object$maximal.mod$x)[[2]])[int_index==1],post_prob=as.numeric(interest),
post_mean=post_mean,post_var=post_var,
lower=lower,upper=upper,thin=thin,prob.level=prob.level,formula=bformula,order=border,
BETA=redBETA,SIG=redSIG,Y0=redY0,TOT=redTOT,meanTOT=meanTOT,int=int,
PRED=PRED,Tpred=Tpred,Tobs=Tobs,pval=pval,statnum=statnum,statistic=statistic)

class(est)<-"submod"

est}
