library(ltm)

itemPar3PLconst<-function(data,c=rep(0,ncol(data))){
J<-ncol(data)
if (length(c)==1) Guess<-rep(c,ncol(data))
else Guess<-c
cont<-cbind(1:ncol(data),rep(1,ncol(data)),Guess)
mod<-tpm(data,constraint=cont,IRT.param=FALSE)
beta<-summary(mod)$coefficients[(J+1):(2*J),1]
alpha<-summary(mod)$coefficients[(2*J+1):(3*J),1]
Sig1<-diag(vcov(mod)[1:J,1:J])
Sig2<-diag(vcov(mod)[(J+1):(2*J),(J+1):(2*J)])
Sig12<-diag(vcov(mod)[1:J,(J+1):(2*J)])
a<-alpha
b<--beta/alpha
c<-Guess
sea<-sqrt(Sig2)
seb<-sqrt(Sig1/alpha^2+beta^2*Sig2/alpha^4-2*beta*Sig12/alpha^3)
seab<-beta*Sig2/alpha^2-Sig12/alpha
par<-cbind(a,b,sea,seb,seab,c)
colnames(par)<-c("a","b","se(a)","se(b)","cov(a,b)","c")
if (is.null(colnames(data))==FALSE) row<-colnames(data)
else{
row<-NULL
for (i in 1:J) row<-c(row,paste("Item",i,sep=""))
}
rownames(par)<-row
return(par)}


