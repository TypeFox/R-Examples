ePCP <-
function(fit,y,alpha=0.05){
fit<-as.matrix(fit)
y<-as.matrix(y)
n<-nrow(fit)
fity<-cbind(fit,y)
out<-matrix(-1,ncol=3)
out[1,1]<-(sum(1-fity[fity[,2]==0,1])+sum(fity[fity[,2]==1,1]))/n
out[1,2]<-out[1,1]-qnorm(1-alpha/2)*sqrt(out[1,1]*(1-out[1,1])/n)
out[1,3]<-out[1,1]+qnorm(1-alpha/2)*sqrt(out[1,1]*(1-out[1,1])/n)
rownames(out)<-c("");colnames(out)<-c("ePCP","lower","upper")
out
}
