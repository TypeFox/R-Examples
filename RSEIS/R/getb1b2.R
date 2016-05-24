`getb1b2` <-
function(J, L, zwin, maxx, max2 )
  {
    wdif = L-J
    b1 = (J-zwin)
    b2 = (L+zwin)
    if(b2>max2) b2 = max2
    if(b1<1) b1 = 1
    
    if((b2-b1)>maxx) { return(c(b1, b2)) }

    while((b2-b1)<maxx)
      {
        b1 = b1-10
        b2 = b2 +10
        if(b2>max2) b2 = max2
        if(b1<1) b1 = 1
        
        if((b2-b1)>maxx) { return(c(b1, b2)) }
 
      }

    return(c(1,1)) 
    
  }

