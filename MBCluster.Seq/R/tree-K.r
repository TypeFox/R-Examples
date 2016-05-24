tree.K=function(k,k1,k2){
  nG=length(k)
  nK=length(k1)
  K=matrix(0,nG,nK)
  K[,1]=k
  for( i in 2:nK){
   k0=K[,i-1]
   k0[k0==k2[i-1]]=k1[i-1]
   K[,i]=k0
 } 
 return(K)
}

