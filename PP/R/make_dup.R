


make_dup <- function(respm)
  {
    
  patts <- apply(respm,1,function(x) paste(x,collapse=""))
  ndpat <- !duplicated(patts)
  seppat <- patts[ndpat]
  posvec <- vector(mode="integer",length = length(patts))
  
  for(i in 1:length(seppat))
    {
      posvec[patts %in% seppat[i]] <- i   
    }
  return(list(posvec=posvec,ndpat=ndpat))
  }  


