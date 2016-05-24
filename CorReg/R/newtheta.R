# ' newtheta for Newton
# ' @param X dataset without intercept
# ' @param Z the structure (square without intercept)
# ' @param B the wheighted structure (p+1)xp
# ' @param Atilde p-sized vector of explicative coeffcients
# ' @param Sigma vector of the standard deviations of the subregressions (p2 long)
# ' @param A vector of the coefficient for the main regression
# ' @param lambda lagrange multiplicators
# ' @param nbit number of iteration for Newton
# ' @export
newtheta<-function(X=X,Z=Z,B=B,Sigma=Sigma,A=A,lambda=NULL,Atilde=Atilde,nbit=1){
   X=cbind(1,X)
   I2=which(colSums(Z)!=0)
   p2=length(I2)
   Z=rbind(0,Z)
   Z[1,I2]=1#on ajoute une constante a chaque ssreg
   Z=cbind(0,Z)
   I2=I2+1
   pz=sum(Z!=0)
   quidroite=which(rowSums(Z[,])!=0)
   p1=length(which(rowSums(Z[-I2,])!=0))#I3 ne doit pas intervenir
   if(is.null(lambda)){
      lambda=rep(1,times=p1)
   }
   n=nrow(X)
   J=matrix(0,ncol=(p1+p2+pz),nrow=(p1+p2+pz))
   barZ=which(Z!=0,arr.ind=T)   
   Fvect=rep(0,times=(p1+p2+pz))   
   for (it in 1:nbit){   
#       print(it)
      for(j in 1:p2){
         I1j=barZ[barZ[,2]==I2[j],1]
         debcolj=nrow(barZ[barZ[,2]<I2[j],])
         colonne=(debcolj+1):(debcolj+sum(Z[,I2[j]])) #sous-reg precedentes+
         Fvect[colonne]=(1/Sigma[j]^2)*t(X[,I1j])%*%(X[,I2[j]]-X[,I1j]%*%B[I1j,I2[j]-1])+A[I2[j]]*lambda[which(Z[quidroite,I2[j]]!=0)]     
         Fvect[pz+p1+j]=Sigma[j]^2-(1/n)*t(X[,I2[j]]-X[,I1j]%*%B[I1j,I2[j]-1])%*%(X[,I2[j]]-X[,I1j]%*%B[I1j,I2[j]-1])
         
         J[pz+p1+j,pz+p1+j]=2*Sigma[j]#bloc J9
         J[pz+p1+j,colonne]=(2/n)*t(X[,I2[j]]-X[,I1j]%*%B[I1j,I2[j]-1])%*%X[,I1j]#bloc J7
         diag(J[pz+(1:p1),which(barZ[,2]==I2[j])])=A[I2[j]]#attention on compte l'intercept #blocJ4
         J[colonne,pz+p1+j]=(-2/(Sigma[j]^3))*t(X[,I1j])%*%(X[,I2[j]]-X[,I1j]%*%B[I1j,I2[j]-1]) #bloc J3 
         diag(J[which(barZ[,2]==I2[j]),pz+(1:p1)])=A[I2[j]]#attention on compte l'intercept #bloc J2
         diag(J[colonne,colonne])=(-1/(Sigma[j]^2))*diag(t(X[,I1j])%*%(X[,I1j]))#bloc J1
      }  
      Fvect[(pz+1):(pz+p1)]=A[quidroite]+as.matrix(B[quidroite,I2-1])%*%A[I2]-Atilde[quidroite]
      if(rcond(J)>10^(-16)){
         matint=solve(J)%*%Fvect
         for (j in 1:p2){
            I1j=barZ[barZ[,2]==I2[j],1]
            debcolj=nrow(barZ[barZ[,2]<I2[j],])
            colonne=(debcolj+1):(debcolj+sum(Z[,I2[j]])) #sous-reg precedentes+
            B[I1j,I2[j]-1]=B[I1j,I2[j]-1]-matint[colonne]
         }
         Sigma=Sigma-matint[-c(1:(pz+p1))]
         lambda=lambda-matint[c((pz+1):(pz+p1))]
      }else{
         print(it)
         print("numerically singular matrix J")
        break
      }    
   }  
   return(list(B=B,Sigma=Sigma,lambda=lambda))
}