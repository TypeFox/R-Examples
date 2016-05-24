cccUst <-
function(dataset,ry,rmet,rtime=NULL,Dmat=NULL,delta=1,cl=0.95){

# Creating data dataset
dades<-data.frame(dataset)
if (length(rtime)==0) dades<-rename.vars(dades,from=c(ry,rmet),to=c("y","met"),info=FALSE) else dades<-rename.vars(dades,from=c(ry,rmet,rtime),to=c("y","met","time"),info=FALSE)


catmet<-unique(dades$met)
if (length(rtime)==0){
valtime<-"1" 
ntime=1}
else {
valtime<-unique(dades$time)
 ntime<-length(valtime)}
ns<-length(dades$y)/(2*ntime)

if (length(Dmat)==0){
Dmat<-diag(rep(1,ntime))}

if(sum(dim(Dmat)==c(ntime,ntime))!=2){
stop("Invalid dimensions in weigth matrix")}

#Sort data
if (length(rtime)>0){
dades<-dades[order(dades$time),]
}

# Creating Y,X data matrix

Y<-array(dades[dades$met== catmet[1],]$y,c(ns,ntime))
X<-array(dades[dades$met== catmet[2],]$y,c(ns,ntime))

colnames(Y)<-valtime
colnames(X)<-valtime


phimat<-array(,c(ns*ns,4))
cont<-0
for (i in 1:ns){
for (j in 1:ns){
cont<-cont+1
phimat[cont,c(3,4)]<-phi(X[i,],Y[i,],X[j,],Y[j,],Dmat,delta) 
phimat[cont,c(1,2)]<-c(i,j)
}
}


colnames(phimat)<-c("i","j","phi1","phi2")

U<-sum((phimat[phimat[,1]!=phimat[,2],3]))/(ns*(ns-1))
V<-sum((phimat[phimat[,1]!=phimat[,2],4]))/(ns*(ns-1))



CCC<-((ns-1)*(V-U))/(U+(ns-1)*V)


phimat1<-phimat[phimat[,1]!=phimat[,2],]
phi1<-tapply(phimat1[,3],phimat1[,1],sum)/(ns-1)
phi2<-tapply(phimat1[,4],phimat1[,1],sum)/(ns-1)

phiv<-cbind(phi1,phi2)
UV<-cbind(U,V)

Saux<-array(0,c(2,2))
C<-array(c(1,0,0,2),c(2,2))
for (i in 1:ns){
Saux<-Saux+t(phiv[i,]-UV)%*%(phiv[i,]-UV)
}

Smat<-C%*%(Saux/(ns^2))%*%C

dev<-array(,c(1,2))
dev[1,1]<-((-1)*ns*(ns-1)*V)/((U+(ns-1)*V)^2)
dev[1,2]<-(ns*(ns-1)*U)/((U+(ns-1)*V)^2)

VarCCC<-dev%*%Smat%*%t(dev)

alpha=1-cl

Z<-0.5*(log((1+CCC)/(1-CCC)))
VarZ<-VarCCC/(((1-CCC)^2)*((1+CCC)^2))
ic.z=Z+c(-1,1)*qnorm(1-alpha/2)*sqrt(VarZ)
ic.ccc=(exp(2*ic.z)-1)/(exp(2*ic.z)+1)
result<-c(CCC,ic.ccc,sqrt(VarCCC),Z,sqrt(VarZ))
conf.lab=paste((1-alpha)*100,"%",sep="")
names(result)<-c("CCC",paste("LL CI",conf.lab),paste("UL CI",conf.lab),"SE CCC","Z","SE Z")
class(result)<-"cccUst"
result
}

