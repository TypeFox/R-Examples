`sepia.colors` <-
function(n, k=1)
  {
    if(missing(k)) k = 1
    
   if(k==1)
     {
       sepia=c(94/255, 38/255, 18/255)
     }
    else
      {
    
    sepia=c(112/255, 66/255, 20/255)
  }
    
    S1 = shade.col(n, sepia, bcol = c(1, 1, 1))
   return(S1)
  }

