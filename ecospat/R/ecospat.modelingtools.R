##
ecospat.maxentvarimport<-function(model,dfvar,nperm){
ref<-predict(model,dfvar)
VarImp<-vector()
for (i in 1:ncol(dfvar)){
print(names(dfvar)[i])
refi<-vector()
for (j in 1:nperm){
cali<-dfvar
cali[,i]<-cali[sample(1:nrow(cali),nrow(cali)),i]
refi<-c(refi,1-cor(ref,predict(model,cali)))
}
VarImp<-c(VarImp,round(mean(refi),3))
}
names(VarImp)<-names(dfvar)
return<-VarImp
}

##

ecospat.Epred<-function(x,w=rep(1,ncol(x)),th=0){

wmean<-function(x,w){
  a<-vector()
  for (i in 1:length(w)){
    a<-c(a,x[i]*w[i])
  }
  return(sum(a)/sum(w))
}

Ebin<-function(vect2bin){
th<-vect2bin[1];vec<-vect2bin[-1]
return(BinaryTransformation(vec,th))
}

if (th!=0){
xi<-rbind(th,x)
xi<-apply(xi,2,Ebin)
Ex<-apply(xi,1,wmean,w)
Ex<-round(cbind(xi,Ex)*1000)
}else{
Ex<-round(apply(x,1,wmean,w))
Ex<-cbind(x,Ex)
}
colnames(Ex)<-c(colnames(x),"E")
return(Ex)
}

dsg<-function(data.coor,mindist){
  xy<-cbind(data.coor,1:nrow(data.coor))
  sample.out<-sample(1:nrow(data.coor),nrow(data.coor))
  z<-xy[sample.out,]
  i<-1
  repeat{
    keep<-z[i,]
    zx<-which(z[,1]>(z[i,1]-mindist)& z[,1]<(z[i,1]+mindist))
    zy<-which(z[,2]>(z[i,2]-mindist)& z[,2]<(z[i,2]+mindist))
    out<-zy[na.exclude(match(zx,zy))]
    if (length(out)>1){
    z<-rbind(keep,z[-out,])}
    #print(i)
    i<-i+1
    if (i>= nrow(z)){break}
  }
return (z[,3])
}

pod<-function(fit,obs,th){
  fitbin<-2*BinaryTransformation(fit,th)
  return(length(which((fitbin+obs)==3))/length(which(obs==1)))
}

specificity<-function(pred,obs){
return(length(which((obs+pred)==0))/(length(which((obs+pred)==0))+length(which((pred-obs)==1))))
}

sensitivity<-function(pred,obs){
return(length(which((obs+pred)==2))/(length(which((obs+pred)==2))+length(which((obs-pred)==1))))
}



ecospat.caleval<-function(data,xy,row.num=1:nrow(data),nrep=1,ratio=0.7,disaggregate=0,pseudoabs=0, npres=0,replace=F){
  SampleMat2<-function (ref, ratio)
  {
    ntot <- length(ref)
    npres <- sum(ref)
    ncal <- ceiling(ntot * ratio)
    pres <- sample(which(ref == 1), ceiling(npres * ratio))
    absc <- sample(which(ref == 0), ncal - length(pres))
    samprows <- list(calibration = c(pres, absc), evaluation = (1:ntot)[-c(pres,
                                                                           absc)])
    return(samprows)
  }
  
  pres<-row.num[sample(which(as.vector(data)==1))]
  abs<-row.num[sample(which(as.vector(data)==0))]
  
  presi<-1:length(pres)
  pres.iter<-c()
  absi<-1:length(abs)
  abs.iter<-c()
  for (i in 1:nrep){
    if (disaggregate !=0){
      pres<-row.num[sample(which(as.vector(data)==1))]
      pres<-pres[dsg(xy,disaggregate)]
    }
    if (npres != 0 & npres < length(pres)){
      if (length(presi)<npres){
        a<-c(presi,sample((1:length(pres))[-presi],npres-length(presi)))
        pres.iter<-cbind(pres.iter,a)
        presi<-c(1:length(pres))[-a]
      }else{a<-sample(presi,npres);pres.iter<-cbind(pres.iter,a);presi<-presi[-a]}
    }else{pres.iter<-cbind(pres.iter,1:length(pres))}
    
    if (pseudoabs != 0 & pseudoabs < length(abs)){
      if (length(absi)<pseudoabs){
        a<- c(absi, sample((1:length(abs))[-absi],pseudoabs-length(absi)))
        abs.iter<-cbind(abs.iter,a)
        absi<-c(1:length(abs))[-a]
      }else{a<-sample(absi,pseudoabs);abs.iter<-cbind(abs.iter,a);absi<-absi[-a]}
    }else{abs.iter<-cbind(abs.iter,1:length(abs))}
    
    y<-rbind(cbind(rep(1,nrow(pres.iter)),pres[pres.iter[,i]]),
             cbind(rep(0,nrow(abs.iter)),abs[abs.iter[,i]]))
    
    cal.eval<-SampleMat2(y[,1],ratio)
    yeval<-y[cal.eval$eval,2]
    ycal<-y[cal.eval$cal,2]
    if (i==1){
      eval<-data.frame(yeval)
      cal<-data.frame(ycal)
    }
    if (i!=1){
      if (length(yeval)>nrow(eval)){eval<-rbind(eval,matrix(nrow=length(yeval)-nrow(eval),ncol=ncol(eval),dimnames=list(c(),rep("yeval",ncol(eval)))))}
      if (length(yeval)<nrow(eval)){yeval<-c(yeval,rep(NA,nrow(eval)-length(yeval)))}
      if (length(ycal)>nrow(cal)){cal<-rbind(cal,matrix(nrow=length(ycal)-nrow(cal),ncol=ncol(cal),dimnames=list(c(),rep("ycal",ncol(cal)))))}
      if (length(ycal)<nrow(cal)){ycal<-c(ycal,rep(NA,nrow(cal)-length(ycal)))}
      eval<-cbind(eval, yeval)
      cal<-cbind(cal,ycal)
    }
  }
  return(list("eval"=eval,"cal"=cal))
}



ecospat.npred<-function(x,th){
x[abs(x)>th]<-NA
count<-function(x){length(which(is.na(x)))}
xi<-x
diag(xi)<-1

while (length(which(is.na(xi))>0)){
narank<-order(apply(xi,1,count),decreasing=T)[1]
xi<-xi[-narank,];xi<-xi[,-narank]
}
return(dim(xi)[1])
}

