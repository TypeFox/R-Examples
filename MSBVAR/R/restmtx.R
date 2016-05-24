"restmtx" <-
function(nsteps,m)
  { tmp <- c(1,rep(0,m-1))
    R <- matrix(0,nrow=nsteps,ncol=(nsteps*m))
    for(i in 1:nsteps)
      { R[i,1:(i*m)] <- matrix(rep(tmp,i),nrow=1) }
    return(R)
  }

