meanTranspose <-
function(x,tol=1e-6)
{
  n = nrow(x)
  p = ncol(x)
  if(sum(is.na(x))==0)
    {
      nu = colMeans(x)
      xnew = scale(x,center=TRUE,scale=FALSE)
      mu = rowMeans(xnew)
      xnew = t(scale(t(xnew),center=TRUE,scale=FALSE))
    }
  else
    {
      ind = 1
      xnew = x
      mut = nut = 0
      while(abs(ind)>tol)
        {
          mu = rowMeans(xnew,na.rm=TRUE)
          xnew = xnew - mu
          nu = colMeans(xnew,na.rm=TRUE)
          xnew = t(t(xnew) - nu)
          ind = sum(mu) + sum(nu)
          mut = mut + mu
          nut = nut + nu
        }
      xnew[is.na(x)] = 0
      mu = mut
      nu = nut
    }
  M = t(t((matrix(0,n,p) + mu)) + nu)
  return(list(x=x,xcen=xnew,mu=mu,nu=nu,M=M))
}

