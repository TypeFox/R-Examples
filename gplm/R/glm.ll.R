"glm.ll" <- function(mu,y,phi=1,family="gaussian", k=1){
  
  "xlny" <- function(x,y){
    x*log(y+(y==0))*(y!=0) + (y==0)*log(1-(x!=0)*(y==0))
  }
  
  if (family=="bernoulli"){
    ll <- xlny(y,mu) + xlny(1-y,1-mu)
  }
  if (family=="binomial"){
    ll <- xlny(y,mu) + xlny(k-y,1-mu)
    ll <- ll + lgamma(k+1)/(lgamma(y+1)*lgamma(k-y+1))
  }
  if (family=="gaussian"){
    ll <- -0.5*(y-mu)^2
    ll <- ll/phi -0.5*log(2*pi) -0.5*log(phi)
  }
  if (family=="poisson"){
    ll <- xlny(y,mu) - mu - lgamma(y+1)
  }
  if (family=="gamma"){
    ll <- -y/mu -log(mu)
    ll <- ll*phi -lgamma(phi) + phi*log(phi*y) -log(y)
  }
  if (family=="inverse.gaussian"){
    ll <- -0.5 * y/(mu*mu) + 1/mu - 0.5/y
    ll <- ll/phi -0.5*log(2*pi) -1.5*log(y) -0.5*log(phi) 
  }
  if (family=="negative.binomial"){
    ll <- y*log(mu) - y*log(mu+k) - k*log(mu+k)
    ll <- ll + lgamma(y+k) - lgamma(y+1) - lgamma(k) + k*log(k)
  }
  return(ll)
}

