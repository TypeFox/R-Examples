"summary.psgc" <-
function(object,...){
qC<-qM.sM(object$C.psamp)
qR<-qM.sM(sR.sC(object$C.psamp))

p<-dim(qC)[1]
nsamp<-dim(object$C.psamp)[3]
vn<-colnames(qC[,,1])
ACF<-QC<-QR<-NULL

rnamesC<-rnamesR<-NULL

for(l1 in 1:p)  {
for(l2 in (1:p)[-l1]) {

QR<-rbind(QR, c(qR[l1,l2,1:3]) )
rnamesR<-c(rnamesR,paste(vn[l1],vn[l2],sep="~")  )

if(l1<l2) {
QC<-rbind(QC, c(qC[l1,l2,1:3]) )
ACF<-rbind(ACF,acf(object$C.psamp[l1,l2,],lag.max=round(nsamp/20),plot=FALSE)[[1]][-1] )
rnamesC<-c(rnamesC, paste(vn[l1],vn[l2],sep="*")  ) 
            }
 }}

rownames(QC)<-rownames(ACF)<-rnamesC
rownames(QR)<-rnamesR

Kappa<-1+2*apply(ACF,1,sum)

ESS<-nsamp/Kappa

ACR.lag1<-matrix(0,p,p)
ACR.lag1[upper.tri(ACR.lag1)]<-ACF[1:choose(p,2),1]
ACR.lag1<-ACR.lag1+t(ACR.lag1)
diag(ACR.lag1)<-NA

nsamp<-dim(object$C.psamp)[3]

SUM<-list(QC=QC,QR=QR,nsamp=nsamp,ACR.lag1=ACR.lag1,ESS=ESS)
structure(SUM,class="sum.psgc")

}

