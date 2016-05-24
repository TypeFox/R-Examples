DE2x = function(y0,times,pars,proc)
{
  parms = list(proc=proc,pars=pars)
  out = lsoda(y0,times,oderhs,parms)
  return(out[,2:ncol(out)])
}


oderhs = function(t,r,parms)
{
  proc = parms$proc
  pars = parms$pars
  r = matrix(r,1,length(r))
  colnames(r) = proc$more$names

  y = as.vector(proc$more$fn(t,r,pars,proc$more$more))

  return(list(y))
}