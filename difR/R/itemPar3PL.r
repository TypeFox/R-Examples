## ESTIMATION OF ITEM RESPONSE MODEL

# require ltm library

library(ltm)

# Items parameters

itemPar3PL<-function(data){
J<-ncol(data)
mod<-tpm(data,IRT.param=FALSE)
gamma<-summary(mod)$coefficients[1:J,1]
beta<-summary(mod)$coefficients[(J+1):(2*J),1]
alpha<-summary(mod)$coefficients[(2*J+1):(3*J),1]
Sig3<-diag(vcov(mod)[1:J,1:J])
Sig1<-diag(vcov(mod)[(J+1):(2*J),(J+1):(2*J)])
Sig2<-diag(vcov(mod)[(2*J+1):(3*J),(2*J+1):(3*J)])
Sig13<-diag(vcov(mod)[1:J,(J+1):(2*J)])
Sig23<-diag(vcov(mod)[1:J,(2*J+1):(3*J)])
Sig12<-diag(vcov(mod)[(2*J+1):(3*J),(J+1):(2*J)])
a<-alpha
b<--beta/alpha
c<-gamma
sea<-sqrt(Sig2)
seb<-sqrt(Sig1/alpha^2+beta^2*Sig2/alpha^4-2*beta*Sig12/alpha^3)
sec<-sqrt(Sig3)
seab<-beta*Sig2/alpha^2-Sig12/alpha
seac<-Sig23
sebc<--Sig13/alpha+beta*Sig23/alpha^2
par<-cbind(a,b,c,sea,seb,sec,seab,seac,sebc)
colnames(par)<-c("a","b","c","se(a)","se(b)","se(c)","cov(a,b)","cov(a,c)","cov(b,c)")
if (is.null(colnames(data))==FALSE) row<-colnames(data)
else{
row<-NULL
for (i in 1:J) row<-c(row,paste("Item",i,sep=""))
}
rownames(par)<-row
return(par)}


