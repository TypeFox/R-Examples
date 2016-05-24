`uberstats` <-
function(box.seq,npts,ninter,d){

  nb <- length(box.seq)

  #ecov = ensemble coverage
  #esup = ensemble support
  #eden = ensemble density
  #npts = total points in dset
  #ninter = total interesting pts in dset
  #intin = total interesting points (#) captured by box sequence
  #totin = total points (#) captured by box sequence


  ecov <- 0
  esup <- 0
  eden <- 0

  for (i in 1:nb){
  
    ecov <- ecov + box.seq[[i]]$marcoverage
    esup <- esup + box.seq[[i]]$marcoverage*ninter/(box.seq[[i]]$y.mean*npts)

  }
  
  eden <- ecov*ninter/(esup*npts)
  
  intin  <- ecov*ninter
  totin  <- esup*npts 

  return(list(ecov=ecov,esup=esup,eden=eden,npts=npts,ninter=ninter,intin=intin,totin=totin,totdims=d))
  
}

