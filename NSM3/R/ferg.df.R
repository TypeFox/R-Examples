ferg.df<-function (x,alpha, mu, npoints,...) 
{
  x = sort(x)
  n = length(x)
  Fn = ecdf(x)
  temp.x = seq(min(x), max(x), length.out=npoints)
  f.n = Fn(temp.x)
  f.0 = mu(temp.x,...)
  
  mu.n=(alpha/(alpha + n))*f.0+(n/(alpha+n))*f.n
  return(mu.n)
}
