# ' Extract best columns of Z in terms of R2
Zorder<-function(Z=Z,X=NULL,orderZ=NULL,p2=1,adj=T,decreasing=T){
   if(is.null(orderZ)){
      R2vect=R2Z(Z=Z,X=X,adj=adj)
      orderZ=order(R2vect,decreasing=decreasing)
   }
   p2=min(sum(colSums(Z)!=0),p2)
   if(p2>0){
      Z[,-orderZ[1:p2]]=0
   }
   return(Z)
}
