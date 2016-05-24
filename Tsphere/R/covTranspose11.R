covTranspose11 <-
function(xc,rhor,rhoc,row=TRUE,sigi.init=NULL,delti.init=NULL,thr=1e-4,maxit=1e3,trace=TRUE,thr.glasso=1e-4,maxit.glasso=1e3,pen.diag=TRUE)
{
  n = nrow(xc)
  p = ncol(xc)
  ndelt = nsig = loglike = NULL
  if(length(delti.init)==0)
    {
      delti = diag(p)
      sigi = diag(n)
    } else{ sigi=sigi.init; delti=delti.init; }
  ind = 1; iter = 0;
  gr = gc = NULL
  if(row){ sr = xc%*%delti%*%t(xc)/p
         }else{
      sc = t(xc)%*%sigi%*%xc/n
      gc = glasso(sc,2*rhoc/n,wi.init=gc$wi,w.init=gc$w,maxit=maxit.glasso,thr=thr.glasso,penalize.diag=pen.diag)
      delti = gc$wi
      ndelt = c(ndelt,(sum(delti!=0) - p)/(p^2 - p))
      sr = xc%*%delti%*%t(xc)/p
    }
  while(ind>thr & iter<maxit)
    {
      iter = iter + 1
      oldS = sigi
      oldD = delti
      gr = glasso(sr,2*rhor/p,wi.init=gr$wi,w.init=gr$w,maxit=maxit.glasso,thr=thr.glasso,penalize.diag=pen.diag)
      sigi = gr$wi
      nsig = c(nsig,(sum(sigi!=0) - n)/(n^2 - n))
      sc = t(xc)%*%sigi%*%xc/n
      gc = glasso(sc,2*rhoc/n,wi.init=gc$wi,w.init=gc$w,maxit=maxit.glasso,thr=thr.glasso,penalize.diag=pen.diag)
      delti = gc$wi
      ndelt = c(ndelt,(sum(delti!=0) - p)/(p^2 - p))
      sr = xc%*%delti%*%t(xc)/p
      ind = (sum((oldD - delti)^2) + sum((oldS - sigi)^2))/(sum((xc^2)))      
      if(trace)
        {
          val = MNloglike(xc,M=matrix(0,n,p),Sig=Sigmahat,Delt=Deltahat,Sigi=sigi,Delti=delti,rhor=rhor,rhoc=rhoc,qr=1,qc=1)
          if(iter>1){if(val<rev(loglike)[1]){ind = 0}}
          loglike = c(loglike,val)
          cat(val,fill=TRUE)
        }
    }
  Sigmahat = gr$w
  Deltahat = gc$w
  return(list(Sigmahat=Sigmahat,Deltahat=Deltahat,Sigmaihat=sigi,Deltaihat=delti,nsig=nsig,ndelt=ndelt,loglike=loglike))
}

