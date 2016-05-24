forward.prediction.error = function(times,data,coefs,lik,proc,pars,whichtimes=NULL)
{
  if(is.null(whichtimes)){ whichtimes = cbind( 1:(length(times)-1),2:length(times) ) }
  
  if(is.matrix(whichtimes)){
    twhich = whichtimes
    whichtimes = list(len=nrow(twhich))
    for(i in 1:nrow(twhich)){
      whichtimes[[i]] = twhich[i,1]:twhich[i,2]
    }
  }

  traj = lik$bvals%*%coefs
  ptraj = c()
  pdata = c()
  ptimes = c()
  pts = c()
  whichobs = c()
  for(i in 1:length(whichtimes)){

      y0 = traj[whichtimes[[i]][1],]
      ts = c(times[whichtimes[[i]]])
  
      parms = list(pars=pars,proc=proc$more)
      out = lsoda(y=y0,times=ts,func=oderhs,parms=parms)

      ptraj = rbind(ptraj,out[2:nrow(out),2:ncol(out)])
      pdata = rbind(pdata,data[whichtimes[[i]][-1],])
      pts  =  c(pts,times[whichtimes[[i]][-1]])
      whichobs = c(whichobs,whichtimes[[i]][-1])
  }

  lik$more$whichobs = whichobs
  return(sum(lik$fn(pdata,pts,ptraj,pars,lik$more)))
}


