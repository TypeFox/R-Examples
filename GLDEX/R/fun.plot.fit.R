"fun.plot.fit" <-
function(fit.obj, data, nclass=50, xlab="", name="", param.vec, ylab="Density", main=""){

fit.obj<-as.matrix(fit.obj)
no.comp<-ncol(fit.obj)
hist.lim<-histsu(data,probability=TRUE,plot=FALSE,xlab="",ylab="",nclass=nclass)$density

jmax<-max(hist.lim,
sapply(1:no.comp,function(i,data,fit.obj,param.vec)
dgl(seq(min(data),max(data),length=1000),fit.obj[1,i],fit.obj[2,i],fit.obj[3,i],fit.obj[4,i],param=param.vec[i]),data,fit.obj,param.vec))

histsu(data,probability=TRUE,xlab=xlab,ylab=ylab,nclass=nclass,ylim=c(0,jmax),main=main)
sapply(1:no.comp,function(i,data,fit.obj,param.vec)
lines(seq(min(data),max(data),length=1000),dgl(seq(min(data),max(data),length=1000),fit.obj[1,i],fit.obj[2,i],fit.obj[3,i],fit.obj[4,i],
param=param.vec[i]),col=i+3,lwd=4),data,fit.obj,param.vec)

name <- rep(name, ncol(fit.obj))
name[dimnames(fit.obj)[[2]] == "STAR"] <- ""

legend("topright",paste(dimnames(fit.obj)[[2]],name,sep=""),lwd=rep(4,no.comp),col=seq(4,3+no.comp))


}

