IntegrateForward = function(y0,ts,pars,proc,more=NULL)
{
  if(is.function(proc))
    proc = list(fn=proc,more=more,names=names(y0))
  else
    proc = proc$more
  parms = list(pars=pars,proc=proc,more=more)
  out = lsoda(y=y0,times=ts,func=oderhs,parms=parms)
  
  return(list(times=out[,1],states=out[,2:ncol(out)]))
}



oderhs = function(t,r,parms)
{

  proc = parms$proc
  pars = parms$pars
 
  r = matrix(r,1,length(r))
  if(!is.null(proc$names))
    colnames(r) = proc$names

  y = as.vector(proc$fn(t,r,pars,proc$more))

  return(list(y))
}