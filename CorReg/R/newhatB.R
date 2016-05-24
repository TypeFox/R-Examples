# ' Computes B based on A and Atilde under the constraint of the subregressions
# ' @param Z the structure p x p binary matrix
# ' @param A the vector of coefficients
# ' @param Atilde the explicative model
# ' @param Bold the old B
# ' @param intercept boolean to define an intercept
# ' @export
newhatB<-function(Z=Z,A=A,Atilde=Atilde,Bold=B,intercept=TRUE){
   qui=WhoIs(Z=Z,I3=F,I2=T,I1=T)
   I2=qui$I2
   B=Bold
   Atilde=matrix(Atilde,ncol=1);A=matrix(A,ncol=1)
   for(j in I2){
      if(as.numeric(A[j+intercept])!=0){#si A2 est nul,  on ne peut utiliser la formule et donc on garde Bold
         I1loc=c(1,which(Z[,j]!=0)+1)
         if(sum(Atilde[I1loc]!=A[I1loc])!=0){#A et Atilde sont egaux donc A2 est nul, on ne peut utiliser la formule et donc on garde Bold
            B[I1loc,j]=(Atilde[I1loc]-A[I1loc])/as.numeric(A[j+intercept])             
            #print("calcul effectif")
         }
      }
   }
   return(B)
}
