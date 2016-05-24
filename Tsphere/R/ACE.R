ACE <-
function(x,sig,delt,sigi,delti,M,thr=1e-4,maxit=1e3)
{
  n = nrow(x)
  p = ncol(x)
  xhat = x
  xhat[is.na(x)] = M[is.na(x)]
  ind = 1
  iter = 1
  rmi = (1:n)[apply(is.na(x),1,sum)>0]
  cmj = (1:p)[apply(is.na(x),2,sum)>0]
  nrmi = apply(is.na(x),1,sum)
  ncmj = apply(is.na(x),2,sum)
  while(ind>thr & iter<maxit)
    {
      oldx = xhat
#by rows
      for(i in rmi)
        {
          ey = M[i,] + (-sigi[i,-i]/sigi[i,i])%*%(xhat[-i,] - M[-i,])
          covy = delt/sigi[i,i]
          mj = (1:p)[is.na(x[i,])]
          if(nrmi[i]>=round(p/2))
            {
              swpz = matinv(covy,(1:p)[-mj])
              xhat[i,mj] = ey[mj] + swpz[mj,-mj,drop=FALSE]%*%(xhat[i,-mj] - ey[-mj])
            }
          else
            {
              swpz = matinv(delti,(1:p)[mj])
              xhat[i,mj] = ey[mj] - swpz[mj,-mj,drop=FALSE]%*%(xhat[i,-mj] - ey[-mj])
            }
        }
#by cols
      for(j in cmj)
        {
          ey = M[,j] + (xhat[,-j] - M[,-j])%*%(-delti[-j,j]/delti[j,j])
          covy = sig/delti[j,j]
          mi = (1:n)[is.na(x[,j])]
          if(ncmj[j]>=round(n/2))
            {
              swpz = matinv(covy,(1:n)[-mi])
              xhat[mi,j] = ey[mi] + swpz[mi,-mi,drop=FALSE]%*%(xhat[-mi,j] - ey[-mi])
            }
          else
            {
              swpz = matinv(sigi*delti[j,j],(1:n)[mi])
              xhat[mi,j] = ey[mi] - swpz[mi,-mi,drop=FALSE]%*%(xhat[-mi,j] - ey[-mi])              
            }
        }
      ind = sum((oldx - xhat)^2)/sum(oldx^2)
      iter = iter + 1
    }
  return(list(x=xhat,iter=iter))
}

