bcct <-
function(formula,data=NULL,n.sample,prior="SBH",start.formula=NULL,start.beta=NULL,start.sig=NULL,save=0,name=NULL,null.move.prob=0.5,a=0.001,b=0.001,progress=FALSE){

if(n.sample<=0){
stop("n.sample must be positive")}
if(prior!="UIP" & prior!="SBH"){
stop("prior not found")}
if(save<0){
stop("save must be non-negative")}
if(null.move.prob<0 | null.move.prob>1){
stop("null.move.prob is a probability and should be between 0 and 1")}
if(a<0 & a!=(-1)){
stop("a and b must be non-negative")}
if(b<0){
stop("a and b must be non-negative")}

ptm<-(proc.time())[3]

if(!is.null(data)){
if(attributes(data)$class=="table"){
data<-data.frame(data)}}

if(save>0){
if(is.null(name)){
name_RJACC<-"RJACC.txt"
name_MHACC<-"MHACC.txt"
name_BETA<-"BETA.txt"
name_MODEL<-"MODEL.txt"
name_SIG<-"SIG.txt"} else{
name_RJACC<-paste(name,"RJACC.txt",sep="")
name_MHACC<-paste(name,"MHACC.txt",sep="")
name_BETA<-paste(name,"BETA.txt",sep="")
name_MODEL<-paste(name,"MODEL.txt",sep="")
name_SIG<-paste(name,"SIG.txt",sep="")}						## names for files if we are saving them

if(file.exists(name_BETA)){stop(paste("A file named ",name_BETA," already exists in the working directory",sep=""))}
if(file.exists(name_MODEL)){stop(paste("A file named ",name_MODEL," already exists in the working directory",sep=""))}
if(file.exists(name_SIG)){stop(paste("A file named ",name_SIG," already exists in the working directory",sep=""))}
if(file.exists(name_RJACC)){stop(paste("A file named ",name_RJACC," already exists in the working directory",sep=""))}
if(file.exists(name_MHACC)){stop(paste("A file named ",name_MHACC," already exists in the working directory",sep=""))}
}										## error if files already exist

priortypes<-c("UIP","SBH")
priornum<-c(1,2)[prior==priortypes]						## which prior?

options(contrasts=c("contr.sum","contr.poly"),warn=-1)
if(!is.null(data)){
maximal.mod<-glm(formula=formula,data=data,family=poisson,control=list(maxit=1),x=TRUE,y=TRUE)} else{
maximal.mod<-glm(formula=formula,family=poisson,control=list(maxit=1),x=TRUE,y=TRUE)}
options(contrasts=c("contr.treatment","contr.poly"),warn=0)			## get maximal design matrix

big.X<-maximal.mod$x								## maximal design matrix
y<-maximal.mod$y									## responses

n<-dim(big.X)[1]								## sample size (number of cells)

IP<-t(big.X)%*%big.X/n
IP[,1]<-0
IP[1,0]<-0									## inverse prior scale matrix	

bmod<-beta_mode(X=big.X,y=y,prior=prior,IP=IP,a=a,b=b)					## find posterior mode

eta.hat<-as.vector(big.X%*%matrix(bmod,ncol=1))					## calculate eta.hat for RJ

if(is.null(start.formula)){
start.index<-rep(1,dim(big.X)[2])} else{
start.index<-formula2index(big.X=big.X,formula=start.formula,data=data)}	## set up starting values for various parameters

if(is.null(start.beta)){
start.beta<-bmod[start.index==1]} 

if(is.null(start.sig)){
start.sig<-1} 

start.mod<-index2model(start.index)

runit<-bcct.fit(priornum=priornum,maximal.mod=maximal.mod,IP=IP,eta.hat=eta.hat,ini.index=start.index,ini.beta=start.beta,ini.sig=start.sig,iters=n.sample,save=save,name=name,null.move.prob=null.move.prob,a=a,b=b,progress=progress)									## do some iterations of data augmentation algorithm
BETA<-runit$BETA
MODEL<-runit$MODEL
SIG<-runit$SIG
rj_acc<-runit$rj_acc
mh_acc<-runit$mh_acc								## list results

if(save>0){
rj_acc<-read.matrix(file=name_RJACC,header=FALSE)
mh_acc<-read.matrix(file=name_MHACC,header=FALSE)
BETA<-read.matrix(file=name_BETA,header=FALSE)
SIG<-read.matrix(file=name_SIG,header=FALSE)
MODEL<-as.character(read.table(file=name_MODEL,header=FALSE)[,1])}

time<-(proc.time())[3]-ptm

est<-list(BETA=BETA,MODEL=MODEL,SIG=SIG,rj_acc=rj_acc,mh_acc=mh_acc,priornum=priornum,maximal.mod=maximal.mod,IP=IP,eta.hat=eta.hat,save=save,name=name,null.move.prob=null.move.prob,time=time,a=a,b=b)

class(est)<-"bcct"								## setup class for results

est}
