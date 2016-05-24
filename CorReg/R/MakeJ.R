#  '  Constructs the jacobian for the constrained likelihood  
#  ' @param B matrice p+1 x p
# ' @param X dataset n x p
# ' @param A Fixed p+1 coefficient vector for the main regression
#  ' @param Z p squared binary matrix describing the structure (adjacency matrix)
# ' @param Sigma p2 vector of standard deviations for each subregression
# '@export
MakeJ<-function(X=X,Z=Z,B=B,Sigma=Sigma,A=A){
   X=cbind(1,X)
   I2=which(colSums(Z)!=0)
   p2=length(I2)
   Z=rbind(0,Z)
   Z[1,I2]=1#on ajoute une constante a chaque ssreg
   Z=cbind(0,Z)
   I2=I2+1
   pz=sum(Z!=0)
   p1=length(which(rowSums(Z[-I2,])!=0))#I3 ne doit pas intervenir
   n=nrow(X)
   J=matrix(0,ncol=(p1+p2+pz),nrow=(p1+p2+pz))
   barZ=which(Z!=0,arr.ind=T)
   for(j in 1:p2){
      I1j=barZ[barZ[,2]==I2[j],1]
      J[pz+p1+j,pz+p1+j]=2*Sigma[j]#bloc J9
      debcolj=nrow(barZ[barZ[,2]<I2[j],])
      colonne=(debcolj+1):(debcolj+sum(Z[,I2[j]])) #sous-reg precedentes+
      J[pz+p1+j,colonne]=(2/n)*t(X[,I2[j]]-X[,I1j]%*%B[I1j,I2[j]-1])%*%X[,I1j]#bloc J7
      diag(J[pz+(1:p1),which(barZ[,2]==I2[j])])=A[I2[j]]#attention on compte l'intercept #blocJ4
      J[colonne,pz+p1+j]=(-2/(Sigma[j]^3))*t(X[,I1j])%*%(X[,I2[j]]-X[,I1j]%*%B[I1j,I2[j]-1]) #bloc J3 
      diag(J[which(barZ[,2]==I2[j]),pz+(1:p1)])=A[I2[j]]#attention on compte l'intercept #bloc J2
      diag(J[colonne,colonne])=(-1/(Sigma[j]^2))*diag(t(X[,I1j])%*%(X[,I1j]))#bloc J1
   }    
   return(J)
}