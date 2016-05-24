ch.ro<-function (x,n,alpha, mu, ...) 
{
  r = length(x)
  k.prime = length(x[x<x[1]])
  l.prime = length(x[x=x[1]])
  g.prime = k.prime + .5*(l.prime+1)
  
  one = ((n+alpha)/(r+alpha))*g.prime
  two = .5*((n-r)/(r+alpha))
  three = (n-r)*((alpha*mu(x[1], ...)))/(r+alpha)
  
  g.hat = one-two+three
  return(g.hat)
  
  
}