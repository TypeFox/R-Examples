CompProbUp <-
function(nAA,nBB,EnAB,prob,MaxHet,vec=NULL) {
  pr <-  prob*4*nAA*nBB/((EnAB+2)*(EnAB+1))
  nvec <- c(vec,pr)
  if(EnAB < MaxHet-2) {
     nvec <- CompProbUp(nAA-1,nBB-1,EnAB+2,pr,MaxHet,nvec)
  }
  return(nvec)
}

