CVcov <-
function(x,maxlam,minlam,steps,pmiss=.01,do=2,trace=TRUE)
{
  n = nrow(x)
  p = ncol(x)
  lams = exp(seq(log(maxlam),log(minlam),l=steps))
  sam = sample(1:(n*p),floor(n*p*pmiss)*do)
  nout = length(sam)/do
  cvmat = matrix(0,do,length(lams))
  for(k in 1:do)
    {
      xx = x
      xx[sam[(nout*(k-1)+1):(nout*k)]] = NA
      mt = meanTranspose(xx)
      xc = mt$xcen
      M = mt$M
      sigi = diag(n)
      delti = diag(p)
      for(i in 1:length(lams))
        {
          ct = covTranspose11(xc,lams[i]*n,lams[i]*p,sigi.init=sigi,delti.init=delti)
          sigi = ct$Sigmaihat
          delti = ct$Deltaihat
          xhat = ACE(xx,ct$Sigmahat,ct$Deltahat,sigi,delti,M)$x
          cvmat[k,i] = mean((xhat[is.na(xx)] - x[is.na(xx)])^2)
          if(trace==TRUE){cat(paste("Fold: ",k,", Lambda: ",lams[i],", Error: ",cvmat[k,i],sep=""),fill=TRUE)}
        }
    }
  optlam = lams[which.min(apply(cvmat,2,mean))]
  return(list(cvmat=cvmat,optlam=optlam,lams=lams))
}

