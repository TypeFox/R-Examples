"fun.plot.many.gld" <-
function(fit.obj,data,xlab="x",ylab="Density",main="",legd="",param.vec){

check.vec<-dimnames(fit.obj)[[2]]

if(length(intersect(check.vec,c("RMFMKL","RPRS")))!=0){
param.vec<-dimnames(fit.obj)[[2]]
param.vec[param.vec=="RPRS"]<-"rs"
param.vec[param.vec!="rs"]<-"fmkl"
no.comp<-ncol(fit.obj)
max.prob<-sapply(1:no.comp, function(i, data, fit.obj, param.vec)
max(dgl(seq(min(data), max(data), length = 1000), fit.obj[1, i], fit.obj[2, i], fit.obj[3, i], fit.obj[4, i], 
param = param.vec[i])),data, fit.obj, param.vec)

matplot(c(min(data),max(data)),c(0,max(max.prob)),xlab=xlab,ylab=ylab,main=main,type="n")
sapply(1:no.comp, function(i, data, fit.obj, param.vec)
lines(seq(min(data), max(data), length = 1000), col=3+i,lwd=4,
dgl(seq(min(data), max(data), length = 1000), fit.obj[1, i], fit.obj[2, i], fit.obj[3, i], fit.obj[4, i], 
param = param.vec[i])),data, fit.obj, param.vec)

legend("topright",legd,lwd=rep(4,no.comp),col=seq(4,3+no.comp))}

if(length(intersect(check.vec,c("RMFMKL","RPRS")))==0){

fit.obj<-matrix(fit.obj,nrow=4)
no.comp <- ncol(fit.obj)
max.prob <- sapply(1:no.comp, function(i, data, fit.obj, 
        param.vec) max(dgl(seq(min(data), max(data), length = 1000), 
        fit.obj[1, i], fit.obj[2, i], fit.obj[3, i], fit.obj[4, 
            i], param = param.vec[i])), data, fit.obj, param.vec)
matplot(c(min(data), max(data)), c(0, max(max.prob)), xlab = xlab, 
        ylab = ylab, main = main, type = "n")
sapply(1:no.comp, function(i, data, fit.obj, param.vec) lines(seq(min(data), 
        max(data), length = 1000), col = 3 + i, lwd = 4, dgl(seq(min(data), 
        max(data), length = 1000), fit.obj[1, i], fit.obj[2, 
        i], fit.obj[3, i], fit.obj[4, i], param = param.vec[i])), 
        data, fit.obj, param.vec)
legend("topright", legd, lwd = rep(4, no.comp), col = seq(4, 
        3 + no.comp))
}

    }