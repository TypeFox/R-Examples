pdsoft.cv <-
function(x, lam.vec=NULL, standard=TRUE, init = c("diag", "soft", "dense"), tau=1e-4, nsplits=10,
                 n.tr=NULL, tolin = 1e-08, tolout = 1e-08, maxitin = 10000, maxitout = 1000, quiet = TRUE)
{
  init=match.arg(init)
  n=dim(x)[1]
  samp.cov=cov(x)*(n-1)/n
  if(is.null(n.tr))
    n.tr=round(n*(1-1/(log(n))))
  n.va=n-n.tr
  
  if(is.null(lam.vec) & standard)
  {
    lam.vec=seq(from=0, to=1, by=0.05)
  }
  if(is.null(lam.vec) & (!standard) )
  {
    samp.minus=abs(samp.cov)
    diag(samp.minus)=0
    lam.vec=seq(from=0, to=max(samp.minus), length.out=20)
  }

  lam.vec=sort(lam.vec, decreasing=TRUE)
  cv.loss = array(0, c(length(lam.vec), nsplits))
  for ( ns in 1:nsplits) 
  {
    ind = sample(n)
    ind.tr = ind[1:n.tr]
    ind.va = ind[(n.tr+1):n]
    s.tr = cov(x[ind.tr,,drop=FALSE]) * (n.tr-1)/n.tr
    s.va = cov(x[ind.va,,drop=FALSE]) * (n.va -1)/n.va
    for( i in 1:length(lam.vec))
    {
      if(i==1)
      {
        out.tmp=pdsoft(s=s.tr, lam=lam.vec[i], tau=tau, init=init, standard=standard, 
                     tolin=tolin, tolout=tolout,
                     maxitin=maxitin, maxitout=maxitout, quiet=quiet)
      } else
      {
        if(standard)
        {
          s0=out.tmp$theta
          i0=out.tmp$theta.inv
        } else
        {
          s0=out.tmp$sigma
          i0=out.tmp$omega
        }
        out.tmp=pdsoft(s=s.tr, lam=lam.vec[i], tau=tau, init="user", s0=s0, 
                     i0=i0,standard=standard, tolin=tolin, tolout=tolout,
                     maxitin=maxitin, maxitout=maxitout, quiet=quiet)
      }
      cv.loss[i,ns] = sum ( (out.tmp$sigma - s.va)^2 )    
      if(!quiet) cat("Finished lam =", lam.vec[i], "in split", ns, "\n") 
    }
    if(!quiet) cat("Finished split", ns, "\n")       
  }  
  cv.err=apply(cv.loss, 1, sum)
  best.lam = lam.vec[which.min(cv.err)]
  out.pds=pdsoft(s=samp.cov, lam=best.lam, tau=tau, init="soft", standard=standard, 
                     tolin=tolin, tolout=tolout,
                     maxitin=maxitin, maxitout=maxitout, quiet=quiet)

  return(list(sigma=out.pds$sigma, omega=out.pds$omega, best.lam=best.lam, cv.err=cv.err, lam.vec=lam.vec, 
              n.tr=n.tr)) 
}

