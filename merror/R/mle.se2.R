"mle.se2" <-
function(v, r, ni)
{
# Computes the squared standard errors for the squared precisions
# before calling this function, compute the MLE's 
# MLE's stored in r
# ni is the no. of items

   N <- length(r)
   
   H <- matrix(0,N,N)
   D <- matrix(0,N,N)
   
   B0 <- vector("numeric",N)
   B1 <- vector("numeric",N)
   
   for(i in 1:N)
   {
     for(j in 1:N)
     {
       if(j!=i) 
       {
         B0[i] <- B0[i] + 1/r[j]
         B1[i] <- B1[i] + v[i,j]/r[j]
       }
     }
   }
   for(i in 1:N)
   {
     H[i,i] <- 0.5*ni*B0[i]^2/(1+B0[i]*r[i])^2  
     for(j in (1:N))
     {
       if(j!=i)
       {
         for(k in 1:N)
         {
           if(k!=i&k!=j) D[i,j] <- D[i,j] + v[j,k]/r[k]
         }
#      cat("\ni =",i," j =",j)
#      print(D)
       H[i,j] <- -(0.5*ni/(r[j]^2*(1+B0[i]*r[i])^2))*
         (1+2*B0[i]*r[i]-((ni-1)/ni)*(B0[i]*v[i,j]+B1[i]-D[i,j]))  
       } 
     }
   }
#  print(B0)
#  print(B1)
#  print(v)
#  print(r)
#  print(ni)
#  print(N)
#  print(H)
   diag(solve(H))
   
}
