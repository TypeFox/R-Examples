fast.band <-
function(x,k)
{
  ### k = number of subdiagonals to keep
  ### x should be centered
  n=dim(x)[1]
  p=dim(x)[2]
  k=min(c(n-2, k, p-1));
  x = as.double(x)
  mode(n) = "integer"
  mode(p) = "integer"
  mode(k) = "integer"
  tmp=matrix(0, nrow=p, ncol=p)
  tmp=as.double(tmp)
  coutput=.C("bchol",xin=x,nin=n, pin=p,kin=k, bcov=tmp)
  band.cov = matrix(coutput$bcov, nrow=p)
  return(band.cov)
}

