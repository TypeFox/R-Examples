diagonalizematrix <-
function(A,ntimes)
{
   nrowA <- nrow(A)
   ncolA <- ncol(A)
  
   Adiag <- matrix(0,nrow=nrowA*ntimes, ncol=ncolA*ntimes)
   
   firsti <- 1
   firstj <- 1
   for (n in 1:ntimes)
   {
      lasti <- firsti+nrowA-1
      lastj <- firstj+ncolA-1
      Adiag[firsti:lasti,firstj:lastj]<-A
      firsti <- lasti+1
      firstj <- lastj+1
   } 
   return (Adiag)
}

