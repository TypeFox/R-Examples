CompProbDown <-
function(nAA,nBB,EnAB,prob,vec=NULL) {
  pr <-  prob*EnAB*(EnAB-1)/(4*(nAA+1)*(nBB+1))
  nvec <- c(vec,pr)
  if(EnAB > 3) {
     nvec <- CompProbDown(nAA+1,nBB+1,EnAB-2,pr,nvec)
  }
  return(nvec)
}

