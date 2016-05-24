bcctu <-
function(object,n.sample,save=NULL,name=NULL,progress=FALSE){

if(n.sample<=0){
stop("n.sample must be positive")}

ptm<-(proc.time())[3]

if(is.null(save)){save<-object$save}
if(is.null(name)){name<-object$name}

if(save>0 & is.null(name)){
name_RJACC<-"RJACC.txt"
name_MHACC<-"MHACC.txt"
name_BETA<-"BETA.txt"
name_MODEL<-"MODEL.txt"
name_SIG<-"SIG.txt"} else{
name_RJACC<-paste(name,"RJACC.txt",sep="")
name_MHACC<-paste(name,"MHACC.txt",sep="")
name_BETA<-paste(name,"BETA.txt",sep="")
name_MODEL<-paste(name,"MODEL.txt",sep="")
name_SIG<-paste(name,"SIG.txt",sep="")}

if(object$save==0 & save>0){

if(file.exists(name_BETA)){stop(paste("A file named ",name_BETA," already exists in the working directory",sep=""))}
if(file.exists(name_MODEL)){stop(paste("A file named ",name_MODEL," already exists in the working directory",sep=""))}
if(file.exists(name_SIG)){stop(paste("A file named ",name_SIG," already exists in the working directory",sep=""))}
if(file.exists(name_RJACC)){stop(paste("A file named ",name_RJACC," already exists in the working directory",sep=""))}
if(file.exists(name_MHACC)){stop(paste("A file named ",name_MHACC," already exists in the working directory",sep=""))}

write.table(file=name_BETA,x=object$BETA,row.names=FALSE,col.names=FALSE,append=TRUE)
write.table(file=name_MODEL,x=object$MODEL,row.names=FALSE,col.names=FALSE,append=TRUE)
write.table(file=name_SIG,x=object$SIG,row.names=FALSE,col.names=FALSE,append=TRUE)
write.table(file=name_RJACC,x=object$rj_acc,row.names=FALSE,col.names=FALSE,append=TRUE)
write.table(file=name_MHACC,x=object$mh_acc,row.names=FALSE,col.names=FALSE,append=TRUE)

}


start.sig<-object$SIG[length(object$SIG)]
start.index<-model2index(object$MODEL[length(object$SIG)],dig=dim(object$BETA)[2])
start.beta<-object$BETA[dim(object$BETA)[1],start.index==1]

runit<-bcct.fit(priornum=object$priornum,maximal.mod=object$maximal.mod,IP=object$IP,eta.hat=object$eta.hat,
ini.index=start.index,ini.beta=start.beta,ini.sig=start.sig,iters=n.sample,save=save,name=name,null.move.prob=object$null.move.prob,
a=object$a,b=object$b,progress=progress)

BETA<-runit$BETA
MODEL<-runit$MODEL
SIG<-runit$SIG
rj_acc<-runit$rj_acc
mh_acc<-runit$mh_acc

if(save>0){
rj_acc<-read.matrix(file=name_RJACC,header=FALSE)
mh_acc<-read.matrix(file=name_MHACC,header=FALSE)
BETA<-read.matrix(file=name_BETA,header=FALSE)
SIG<-read.matrix(file=name_SIG,header=FALSE)
MODEL<-as.character(read.table(file=name_MODEL,header=FALSE)[,1])}

if(save==0){
rj_acc<-c(object$rj_acc,rj_acc)
mh_acc<-c(object$mh_acc,mh_acc)
BETA<-rbind(object$BETA,BETA)
SIG<-c(object$SIG,SIG)
MODEL<-c(object$MODEL,MODEL)}

ptm<-(proc.time())[3]-ptm
time<-object$time+ptm

est<-list(BETA=BETA,MODEL=MODEL,SIG=SIG,rj_acc=rj_acc,mh_acc=mh_acc,priornum=object$priornum,maximal.mod=object$maximal.mod,IP=object$IP,eta.hat=object$eta.hat,save=save,name=name,null.move.prob=object$null.move.prob,time=time,a=object$a,b=object$b)

class(est)<-"bcct"

est}
