
pred.response = function(fit, newx=NULL)
{
  if(is.null(newx))
  {
    R=t(fit$Xc%*%fit$Rmat)
  } else
  {
    newxc=scale(newx, center=fit$mx, scale=FALSE)
    R=t(newxc%*%fit$Rmat)
  }
  GWG=fit$GWG
  y=fit$y
  BFt=tcrossprod(fit$bhat, fit$fc)
  n=length(y)
  d=dim(GWG)[1]
  nx=dim(R)[2]
  Ehat=rep(0,nx)
  ghat=rep(0, n)
  coutput=.C("pred",  Ehat=as.double(Ehat), ghat=as.double(ghat), R=as.double(R), 
  BFt=as.double(BFt),  GWG=as.double(GWG), y=as.double(y), nin=as.integer(n), din=as.integer(d),
  nxin=as.integer(nx))
  return(coutput$Ehat)
}

  