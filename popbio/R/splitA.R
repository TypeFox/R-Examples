splitA<-function(A, r=1, c=-1)
{
   tm<-A
   fm<-A
   ### check if r is a matrix
   if(is.matrix(r) ){
     # logical matrix (TRUE=fertility)
     if(is.logical(r)){
       tm[r]  <- 0
       fm[!r] <- 0    
     } else{
      tm <- A - r
      fm <- r
     }
   } else{
     tm[r, c] <- 0
     fm[-r, ] <- 0
     fm[r, -(c)] <- 0
   }   
   list(T=tm,F=fm)
}

