######################################

Hybrid.Tree.Microarray=function(data,cluster0,distance="euclidean"){
   X=data$logFC
   k=cluster0
   nK=length(unique(cluster0))
   nG=nrow(X)
   k1=k2=d=rep(0,nK-1)
   for(i in 1:(nK-1)){
     print(paste("level",i))
     D=dst2center.pairs(X,k,distance)
     k1[i]=D[which.min(D[,3]),1]
     k2[i]=D[which.min(D[,3]),2]
     d[i]=min(D[,3])
     k[k==k2[i]]=k1[i]
   } 
     return(cbind(k1,k2,d))
}




dst2center.pairs=function(X,k,distance){
 k1=k2=unique(k)
 nK=length(k1)
 if(nK==1) return(c(k1,k2,0))
 D=c()
 for(i in 1:(nK-1))
  for(j in (i+1):nK){
   x1=X[k==k1[i],]
   x2=X[k==k2[j],]
   if(distance=="euclidean") d=dst.euclidean.pair(x1,x2)
   if(distance=="pearson") d=dst.pearson.pair(x1,x2)
   if(distance=="maximum") d=dst.maximum.pair(x1,x2)
   D=rbind(D,c(k1[i],k2[j],d))
  }
 D=matrix(D,ncol=3)
 dimnames(D)[[2]]=c("k1","k2","D")
 return(D)
}


dst.euclidean.pair=function(x1,x2){
  x0=rbind(x1,x2)
  x1=matrix(x1,ncol=ncol(x0))
  x2=matrix(x2,ncol=ncol(x0))
  x0=sweep(x0,2,meanCol(x0))
  x1=sweep(x1,2,meanCol(x1))
  x2=sweep(x2,2,meanCol(x2))
  d=sqrt(sum(x0^2)-sum(x1^2)-sum(x2^2))
  return(d)
}

dst.pearson.pair=function(x1,x2){
  x0=rbind(x1,x2)
  x1=matrix(x1,ncol=ncol(x0))
  x2=matrix(x2,ncol=ncol(x0))
  x0=x0-meanRow(x0)
  x1=x1-meanRow(x1)
  x2=x2-meanRow(x2)
  z1=x1/sqrt(sumRow(x1^2))
  z2=x2/sqrt(sumRow(x2^2))
  z0=x0/sqrt(sumRow(x0^2))
  m1=meanCol(z1)
  m2=meanCol(z2)
  m0=meanCol(z0)
  m1=m1/sqrt(sum(m1^2))
  m2=m2/sqrt(sum(m2^2))
  m0=m0/sqrt(sum(m0^2))
  d=sum(1-z0%*%m0)-sum(1-z1%*%m1)-sum(1-z2%*%m2)
  return(d)
}
dst.maximum.pair=function(x1,x2){
  x0=rbind(x1,x2)
  x1=matrix(x1,ncol=ncol(x0))
  x2=matrix(x2,ncol=ncol(x0))
  x0=sweep(x0,2,meanCol(x0))
  x1=sweep(x1,2,meanCol(x1))
  x2=sweep(x2,2,meanCol(x2))
  d=sum(maxRow(x0))-sum(maxRow(x1))-sum(maxRow(x2))
  return(d)
}


