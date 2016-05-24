library(ltm)

itemPar2PL<-function(data){
J<-ncol(data)
mod<-ltm(data~z1,IRT.param=FALSE)
beta<-summary(mod)$coefficients[1:J,1]
alpha<-summary(mod)$coefficients[(J+1):(2*J),1]
Sig1<-diag(vcov(mod))[1:J]
Sig2<-diag(vcov(mod))[(J+1):(2*J)]
Sig12<-diag(vcov(mod)[1:J,(J+1):(2*J)])
a<-alpha
b<--beta/alpha
sea<-sqrt(Sig2)
seb<-sqrt(Sig1/alpha^2+Sig2*beta^2/alpha^4-2*beta*Sig12/alpha^3)
sab<-beta*Sig2/alpha^2-Sig12/alpha
par<-cbind(a,b,sea,seb,sab)
colnames(par)<-c("a","b","se(a)","se(b)","cov(a,b)")
if (is.null(colnames(data))==FALSE) 
row<-colnames(data)
else {
row<-NULL
for (i in 1:J) row<-c(row,paste("Item",i,sep=""))
}
rownames(par)<-row
return(par)
}