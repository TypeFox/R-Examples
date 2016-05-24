
FGalgorithm<-function(eF,eG,p,n,A){
  
  Phi<-function(B,A,n){#A is a list of positive definite symetric (p.d.s) matrixs
    K<-length(A);likelihood<-numeric(K)
    for(i in 1:K){likelihood[i]<-((det(diag(diag(t(B)%*%A[[i]]%*%B),length(diag(t(B)%*%A[[i]]%*%B))))/(det(t(B)%*%A[[i]]%*%B)))^(n[i]))}
    return(prod(likelihood))}
    FGPhi<-function(B,A,n, type = 'local'){
    K<-length(A);likelihood<-numeric(K)
    for(i in 1:K){
    rr<-t(B)%*%A[[i]]%*%B
    likelihood[i]<-((det(diag(diag(rr),length(diag(rr))))/(det(rr)))^(n[i]))}
    return(prod(likelihood))
  }
    
 Galgo<-function(eG,n,Tchain, type = 'local'){
    Qchain<-vector("list");K<-length(Tchain)
    Delta<-matrix(nrow=2,ncol=K);T<-vector("list",length=K)
    Q<-diag(2);g<-1
    repeat{
      Qchain[[g]]<-Q
      for(j in 1:2){for(i in 1:K){Delta[j,i]<-t(Q[,j])%*% Tchain[[i]]%*%Q[,j]}}
      for(i in 1:K){
        T[[i]]<-(n[i])*((Delta[1,i]-Delta[2,i])/(Delta[1,i]*Delta[2,i]))*Tchain[[i]]}
      S=mat.or.vec(2, 2);for(i in 1:K){S<-(S+T[[i]])}
      Q[,1]<-(eigen(S)$vectors)[,which(eigen(S)$values==max(eigen(S)$values))]
      if(Q[1,1]<0){Q[,1]<-(-Q[,1])}
      Q[,2]<-(eigen(S)$vectors)[,which(eigen(S)$values==min(eigen(S)$values))]
      if(Q[1,2]<0){Q[,2]<-(-Q[,2])}
      if(norm(Qchain[[g]]-Q,"F")<=eG){break}
      g<-g+1
    }
    return(Q)
  }

  K<-length(A);Tchain<-vector("list",length=K);B<-diag(p)
  Bchain<-vector("list")
  f<-1
  repeat{
    Bchain[[f]]<-B
   for(l in 1:(p-1)){
     for(j in (l+1):p){
        for(i in 1:K){Tchain[[i]]<-(t(B[,c(l,j)])%*%A[[i]]%*%B[,c(l,j)])}
        Q<-Galgo(eG,n,Tchain)
        C<-B[,c(l,j)]%*%Q[,1];D<-B[,c(l,j)]%*%Q[,2]
        B[,l]<-C;B[,j]<-D
      }
    }
    gozare<-(abs((Phi(Bchain[[f]],A,n)-Phi(B,A,n))<eF))
    if(gozare){break}
    f<-f+1
  }
  return(B)
}