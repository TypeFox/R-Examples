ptest=function(X,Xc, y, f, Bhat, weight.matrix, fc, Phi, W, d0, NPERM=1e3, qrtol)
{
  n=nrow(X)
  p=ncol(X)
  if(d0 > 0)  ## used to get Gamhat under null, which is nill when d0=0
  { 
    out.null=solve.exact(Phi=Phi, Bhat=Bhat, d=d0, W=W, qrtol=qrtol)
    Pgw=out.null$G%*%t(out.null$Rmat)
    Qgw=diag(p)-Pgw
  }
  ## get test statistic
  out.alt=solve.exact(Phi=Phi, Bhat=Bhat, d=(d0+1), W=W, qrtol=qrtol)
  R.alt  = Xc%*%out.alt$Rmat
  yhat=pred.pfc(R=t(R.alt), BFt=tcrossprod(out.alt$B, fc), GWG=out.alt$GWG, y=y)
  test.stat=mean((y-yhat)^2)
  
  nulldist=NULL
  for(i in 1:NPERM)
  {
    if( d0 > 0)
    {
      xperm=X%*%t(Pgw) + (X%*%t(Qgw))[sample(n),]
    }else
    {
      xperm=X[sample(n),]
    }  
    permfit=fit.pfc(X=xperm, y=y, d=(d0+1), F.user=f, weight.type=weight.matrix, qrtol=qrtol)
    nulldist[i]=mean((y-pred.response(permfit))^2)
  }
  return(list(stat=test.stat, ndist=nulldist))
}
