EPV<-function (itemBank, item, x, theta, it.given, model=NULL, priorDist = "norm", 
    priorPar = c(0, 1), D = 1, parInt = c(-4,4, 33)) {
th <- theta
if (is.null(model)){
itj <- rbind(it.given, itemBank[item,])
p1 <- Pi(th,itemBank[item,], D = D)$Pi
p0 <- 1 - p1
th0<-eapEst(itj,c(x,0),D=D,priorDist=priorDist,priorPar=priorPar,lower=parInt[1],upper=parInt[2],nqp=parInt[3])
th1<-eapEst(itj,c(x,1),D=D,priorDist=priorDist,priorPar=priorPar,lower=parInt[1],upper=parInt[2],nqp=parInt[3])
var0<-(eapSem(th0,itj,c(x,0),D=D,priorDist=priorDist,priorPar=priorPar,lower=parInt[1],upper=parInt[2],nqp=parInt[3]))^2
var1<-(eapSem(th1,itj,c(x,1),D=D,priorDist=priorDist,priorPar=priorPar,lower=parInt[1],upper=parInt[2],nqp=parInt[3]))^2
res<-p0*var0+p1*var1
}
else{
probs<-Pi(th,itemBank[item,],model=model,D=D)$Pi
probs<-probs[!is.na(probs)]
it.new<-rbind(it.given,itemBank[item,])
th.new<-se.new<-NULL
for (j in 0:(length(probs)-1)){
th.new[j+1]<-eapEst(it=it.new,x=c(x,j),model=model,D=D,priorDist=priorDist,priorPar=priorPar,lower=parInt[1],upper=parInt[2],nqp=parInt[3])
se.new[j+1]<-(eapSem(th.new[j+1],it=it.new,x=c(x,j),D=D,model=model,priorDist=priorDist,priorPar=priorPar,lower=parInt[1],upper=parInt[2],nqp=parInt[3]))^2
}
res<-sum(probs*se.new)
}
return(as.numeric(res))}