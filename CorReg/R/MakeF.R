# Atilde de taille p1
# ' @export
MakeF<-function(X=X,Z=Z,B=B,Sigma=Sigma,A=A,lambda=NULL,Atilde=Atilde){
   X=cbind(1,X)
   I2=which(colSums(Z)!=0)
   p2=length(I2)  
   Z=rbind(0,Z)
   Z[1,I2]=1#on ajoute une constante a chaque ssreg
   Z=cbind(0,Z)
   I2=I2+1
   pz=sum(Z!=0)
   quidroite=which(rowSums(Z[,])!=0)#pour lambda
   p1=length(which(rowSums(Z[-I2,])!=0))#I3 ne doit pas intervenir
   if(is.null(lambda)){
      lambda=rep(1,times=p1)
   }
   n=nrow(X)
   Fvect=rep(0,times=(p1+p2+pz))
   barZ=which(Z!=0,arr.ind=T)
   for(j in 1:p2){
      I1j=barZ[barZ[,2]==I2[j],1]
      debcolj=nrow(barZ[barZ[,2]<I2[j],])
      colonne=(debcolj+1):(debcolj+sum(Z[,I2[j]])) #sous-reg precedentes+
      Fvect[colonne]=(1/Sigma[j]^2)*t(X[,I1j])%*%(X[,I2[j]]-X[,I1j]%*%B[I1j,I2[j]-1])+A[I2[j]]*lambda[which(Z[quidroite,I2[j]]!=0)]     
      Fvect[pz+p1+j]=Sigma[j]^2-(1/n)*t(X[,I2[j]]-X[,I1j]%*%B[I1j,I2[j]-1])%*%(X[,I2[j]]-X[,I1j]%*%B[I1j,I2[j]-1])    
   }  
   Fvect[(pz+1):(pz+p1)]=A[quidroite]+as.matrix(B[quidroite,I2-1])%*%A[I2]-Atilde[quidroite]
   return(Fvect)
}