
MEI<-function (itemBank, item, x, theta, it.given, model=NULL, method="BM", 
     priorDist = "norm", priorPar = c(0, 1),D=1, range=c(-4,4),parInt=c(-4,4,33), infoType="observed") 
{
if (infoType!="Fisher" & infoType!="observed") stop("'infoType' must be either 'Fisher' or 'observed'",call.=FALSE)
   th<-theta
if (is.null(model)){
   itj<-rbind(it.given,itemBank[item,])
   th0<-thetaEst(itj,c(x,0),D=D,method=method,priorDist=priorDist,priorPar=priorPar,range=range,parInt=parInt)
   th1<-thetaEst(itj,c(x,1),D=D,method=method,priorDist=priorDist,priorPar=priorPar,range=range,parInt=parInt)
   p1<-Pi(th,itemBank[item,],D=D)$Pi
   p0<-1-p1
if (infoType=="Fisher"){
Ij0<-Ii(th0,rbind(itemBank[item,]),D=D)$Ii
Ij1<-Ii(th1,rbind(itemBank[item,]),D=D)$Ii
}
else{
   Ij0<-OIi(th0,rbind(itemBank[item,]),0,D=D)
   Ij1<-OIi(th1,rbind(itemBank[item,]),1,D=D)
}
res<-p0*Ij0+p1*Ij1
}
else{
probs<-Pi(th,itemBank[item,],model=model,D=D)$Pi
probs<-probs[!is.na(probs)]
it.new<-rbind(it.given,itemBank[item,])
I.new<-NULL
for (j in 0:(length(probs)-1)){
th.new<-thetaEst(it=it.new,x=c(x,j),model=model,method=method,priorDist=priorDist,priorPar=priorPar,range=range,parInt=parInt,D=D)
if (infoType=="Fisher") I.new<-c(I.new,Ii(th.new,itemBank[item,],model=model,D=D)$Ii)
else I.new<-c(I.new,OIi(th.new,itemBank[item,],x=j,model=model,D=D))
}
res<-sum(probs*I.new)
}
return(as.numeric(res))}