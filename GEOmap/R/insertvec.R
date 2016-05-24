insertvec<-function(v, ind, val)
  {
    sects = list()
    isects = length(ind)+1
    JIND  = c(0, ind, length(v))
    
    for(k in 1:isects)
      {
        J1 = JIND[k]+1
        J2 = JIND[k+1]
        
    sects[[k]] = v[J1:J2]
      }

    w = sects[[1]]
    for(k in 2:(isects))
      {
        w = c(w, val, sects[[k]])
      }
    

    
    return(w)
  }

