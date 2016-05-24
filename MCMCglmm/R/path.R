"path"<-function(cause=NULL, effect=NULL, k){
 
   if(length(cause)!=length(effect)){stop("cause and effect must be the same length")}
   if(all(cause%%1!=0)){stop("cause must be integers")}
   if(all(effect%%1!=0)){stop("effect must be integers")}

   psi<-matrix(0, k, k)
   for(i in 1:length(effect)){
     psi[effect[i],cause[i]]<-1
   }
   psi
}


