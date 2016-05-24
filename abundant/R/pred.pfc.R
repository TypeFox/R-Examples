
pred.pfc = function(R,BFt, GWG, y)
{
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
