pomp.skeleton = function(times,y,p,more)
{
  # Turns a skeleton function from a `pomp' object into the right hand side
  # of an ODE for use in CollocInfer

#  x = array(t(y),dim=c(ncol(y),1,nrow(y)),dimnames=list(colnames(y),NULL,NULL))
#  params = array(data=p,dim=c(length(p),1),dimnames=list(names(p),NULL))

#  y = skeleton(more$pompo.obj,x=y,params=params,t=times)
#  return(t(y[,1,]))
}

pomp.dmeasure = function(times,data,x,p,more)
{
#  x = array(t(x),dim=c(ncol(x),1,nrow(x)),dimnames=list(colnames(x),NULL,NULL))
#  params = array(data=p,dim=c(length(p),1),dimnames=list(names(p),NULL))
#  f = dmeasure(more$pomp.obj,y=t(data),x=x,params=params,times=times,log=TRUE)
#  return(t(f))
}

