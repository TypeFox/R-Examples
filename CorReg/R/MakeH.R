# 'B matrice p+1 x p   

MakeH<-function(X=X,Z=Z,B=B,Sigma=Sigma){
   X=cbind(1,X)
   I2=which(colSums(Z)!=0)
   Z=rbind(0,Z)
   Z[1,I2]=1#on ajoute une constante a chaque ssreg
   Z=cbind(0,Z)
   I2=I2+1
   p2=length(I2)
   pz=sum(Z!=0)
   p1=ncol(X)-p2#prendre donc en compte la constante
   n=nrow(X)
   H=matrix(0,ncol=(p1+p2+pz),nrow=(p1+p2+pz))
   barZ=which(Z!=0,arr.ind=T)
   for(j in 1:p2){
      I1j=barZ[barZ[,2]==I2[j],1]
      H[j,j]=-n/(Sigma[j]^2)+(1/(Sigma[j]^4)) * t(X[,I2[j]]-X[,I1j]%*%B[I1j,I2[j]])%*%(X[,I2[j]]-X[,I1j]%*%B[I1j,I2[j]])   
      debcolj=nrow(barZ[barZ[,2]<I2[j],])
      colonne=(debcolj+1):(debcolj+sum(Z[,I2[j]])) #sous-reg precedentes+
      colonne=colonne+p2#on decale du bloc H1
      H[j,colonne]=(-2/Sigma[j]^3)*t(X[,I1j])%*%(X[,I2[j]]-X[,I1j]%*%B[I1j,I2[j]])
      H[colonne,j]= H[j,colonne]
      H[colonne,colonne]=(-1/(Sigma[j]^2)) * t(X[,I1j])%*%X[,I1j]
   }
   
#    require(matrixcalc)
#    print(is.negative.definite(H))
#    print(is.positive.definite(H))
   return(H)
}