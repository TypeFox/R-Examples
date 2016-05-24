rerow <-
function(w){

 # Here we normalize each row
 
 m=nrow(w)
 for(i in 1:m){
  w[i,]=w[i,]/(sqrt(sum(w[i,]^2)))
 }
 perm = 1:m;

 # Here we reorder the rows
     
 for (j in 1:(m-1)){
  mx=abs(w[j,j])
  kmx=j
  for( k in j:m){
   if (abs(w[k,j])>mx){
    kmx=k
    mx=abs(w[k,j])
   }
  }

  if (kmx!=j){
   y=w[j,]
   w[j,]=w[kmx,]
   w[kmx,]=y 
   t_perm = perm[j]
   perm[j] = perm[kmx]
   perm[kmx] = t_perm        
  }

  if ( w[j,j]<0){
   w[j,]=-w[j,]
  }
 }
     
 if (w[m,m]<0){
  w[m,]=-w[m,]
 }

 nw=w
 return(nw)

}
