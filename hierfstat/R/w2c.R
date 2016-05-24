#########################################################################
wc2<-function(ndat,diploid=TRUE,pol=0.0){
#specific for 2 populations only
cl<-match.call()  
if (!diploid) {
dum<-ndat[,-1]
nd<-max(dum,na.rm=T)
modu<-1000
if (nd<10) modu<-10
if(nd<100) modu<-100
dum<-dum*modu+dum
ndat<-data.frame(ndat[,1],dum)
}

pop<-ndat[,1]
ni<-length(pop)
dat<-ndat
loc.names<-names(dat)[-1]
n<-t(ind.count(dat))
nt<-apply(n,1,sum,na.rm=TRUE)
untyped.loc<-which(nt==0)
typed.loc<-which(nt!=0)

if (length(untyped.loc)>0){
  dat<-dat[,-(untyped.loc+1)]
  n<-t(ind.count(dat))
  nt<-apply(n,1,sum,na.rm=TRUE)
}

alploc<-nb.alleles(cbind(rep(1,ni),dat[,-1]))
np<-2
npl<-apply(n,1,function(x) sum(!is.na(x))) #accounts for pops not typed for one locus
nl<-dim(n)[1]
p<-pop.freq(dat,diploid)
pb<-pop.freq(cbind(rep(1,ni),dat[,-1]),diploid)
n<-matrix(unlist(n),ncol=np) #why?
nal<-n[rep(1:nl,alploc),]
nc<-(nt-apply(n^2,1,sum,na.rm=TRUE)/nt)/(npl-1) #what happens if npl=1?
ntal<-rep(nt,alploc)
ncal<-rep(nc,alploc)
p<-matrix(unlist(lapply(p,t)),ncol=np,byrow=TRUE)
pb<-matrix(unlist(pb),ncol=1)


if (diploid){
dum<-getal.b(dat[,-1])
all.loc<-apply(dum,2,tempfun1<-function(y) as.numeric(dimnames(table(y))[[1]]))

hetpl<-apply(dum,2,fun<-function(z){
lapply(as.numeric(dimnames(table(z))[[1]]),
who.is.het<-function(y) apply(z==y,1,
ind.is.het<-function(x) xor(x[1],x[2])))}
)

mho<-lapply(hetpl,
tempfun2<-function(x) matrix(unlist(lapply(x,
tempfun3<-function(y) tapply(y,pop,sum,na.rm=TRUE))),ncol=np)
)
mho<-matrix(unlist(mho),ncol=np,byrow=TRUE)
mhom<-(2*nal*p-mho)/2
}
else mhom<-nal*p

SSG<-apply(nal*p-mhom,1,sum,na.rm=TRUE)

dum<-nal*(p-2*p^2)+mhom
SSi<-apply(dum,1,sum,na.rm=TRUE)

dum1<-nal*(sweep(p,1,pb))^2
SSP<-2*apply(dum1,1,sum,na.rm=TRUE)

ntalb<-rep(npl,alploc)   #corrects for untyped samples at a locus

MSG<-SSG/ntal
MSP<-SSP/(ntalb-1)
MSI<-SSi/(ntal-ntalb)

sigw<-MSG
sigb<-0.5*(MSI-MSG)
siga<-1/2/ncal*(MSP-MSI)

Fxy<-function(x) x[1]/sum(x,na.rm=TRUE)

FST.pal<-apply(cbind(siga,sigb,sigw),1,Fxy)
FIS.pal<-apply(cbind(sigb,sigw),1,Fxy)

loc<-rep(1:nl,alploc)
lsiga<-tapply(siga,loc,sum,na.rm=TRUE)
lsigb<-tapply(sigb,loc,sum,na.rm=TRUE)
lsigw<-tapply(sigw,loc,sum,na.rm=TRUE)

sigloc<-cbind(lsiga,lsigb,lsigw)

if (length(untyped.loc)>0){
x<-order(c(typed.loc,untyped.loc))
dum1<-matrix(numeric(3*length(untyped.loc)),ncol=3)
dum<-rbind(sigloc,dum1,deparse.level=0)[x,]

sigloc<-dum
}
sigloc<-data.frame(sigloc)
names(sigloc)=c("lsiga","lsigb","lsigw")
rownames(sigloc)<-loc.names

lFST<-apply(cbind(lsiga,lsigb,lsigw),1,Fxy)
lFIS<-apply(cbind(lsigb,lsigw),1,Fxy)

tsiga<-sum(siga,na.rm=TRUE)
tsigb<-sum(sigb,na.rm=TRUE)
tsigw<-sum(sigw,na.rm=TRUE)
tFST<-Fxy(c(tsiga,tsigb,tsigw))
tFIS<-Fxy(c(tsigb,tsigw))

res<-list(call=cl,sigma=cbind(loc,siga,sigb,sigw),
sigma.loc=sigloc,
per.al=list(FST=FST.pal,FIS=FIS.pal),
per.loc=list(FST=lFST,FIS=lFIS),
FST=tFST,FIS=tFIS)

if (!diploid){
res<-list(call=cl,sigma=cbind(loc,siga,sigb,sigw),
sigma.loc=sigloc,
per.al=list(FST=FST.pal),
per.loc=list(FST=lFST),
FST=tFST,FIS=NA)
}
class(res)<-"wc"
res
}

print.wc<-function(x,...){
print(list(FST=x$FST,FIS=x$FIS))
invisible(x)
}

