khmat = function(levsind,permh){

# jk 06-26-2011 make one-line functions internal
  avea<-function(a) matrix(1/a,ncol=a,nrow=1)

  deltaa = function(a){
     am1 = a - 1 ; cbind(diag(am1),rep(-1,am1))
  }


   k = length(levsind[1,])
   n = length(levsind[,1])

   numba = length(unique(levsind[,1]))
   if(permh[1] == 1){
       khmat = deltaa(numba)
   } else {
       khmat = avea(numba)
   }
   for(i in 2:k){
      numba = length(unique(levsind[,i]))
      if(permh[i] == 1){
       khmat = kronecker(khmat,deltaa(numba))
      } else {
          khmat = kronecker(khmat,avea(numba))
      }
   }
   khmat
}
