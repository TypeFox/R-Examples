somme<-function(X,poids)
{
  return(sum(X*poids))
}


contrindiv<-function(Pinoyau)
{
  ncrit<-ncol(Pinoyau)
  nprod<-nrow(Pinoyau)
  vpi<-vector("list")
  u<-vector("list")
  denom<-vector("list")
  res<-vector("list")
  for (k in 1:ncrit)
  {
    vpi[[k]]<-rbind(t(Pinoyau[,k]),matrix(Pinoyau[,k],nrow=nprod,ncol=nprod))
    u[[k]]<-diag(1,nprod)
    u[[k]]<-cbind(rep(1,nprod),u[[k]])
    denom[[k]]<-u[[k]]%*%vpi[[k]]
    vpi[[k]]<-t(matrix(Pinoyau[,k],nrow=nprod,ncol=nprod))
    res[[k]]<-t(vpi[[k]]/denom[[k]])
  }
  return(res)
}


estimlambda<-function(Pi,Data,lambda)
{
  Tcla<-length(Pi)
  ncrit<-length(Data@Paircomp)
  nsujet<-length(Data@Cons)
  puissance<-vector("list")
  for (t in 1:Tcla)
  {
    contrindivt<-contrindiv(Pi[[t]])
    puissancet<-vector("list")
    for (k in 1:ncrit)
    {
      base<-as.matrix(contrindivt[[k]])
      Matbase<-lapply(Data@Paircomp[[k]],as.matrix)
      puissancet[[k]]<-simplify2array(lapply(lapply(Matbase,contri,base),prod))
    }
    puissance[[t]]<-exp(rowSums(log(simplify2array(puissancet))))
    puissance[[t]]<-puissance[[t]]*lambda[t]
  }
  lvrnouv<-sum(log(rowSums(simplify2array(puissance))))
  Zhtnew<-(simplify2array(puissance))/(rowSums(simplify2array(puissance)))
  lambdanew<-colSums(Zhtnew)/nsujet
  resu<-list(lambdanouv=lambdanew,Zhtnouv=Zhtnew,lvrnouv=lvrnouv)
  return(resu)
}


contri<-function(X,A)
{
  return(A^X)
}
