nfix <-
function(r, l, u, theta, sigma, type1, type2, nrange=c(0,1000))
{ 
  f <- function(x, theta, sigma, type1, type2, l,u){
    sigma.of.est.theta <- sigma*sqrt(1/x +1/(r*x))
    pnorm((u - theta)/(sigma.of.est.theta)- qnorm(1-type1))- 
    pnorm((l - theta)/(sigma.of.est.theta)+ qnorm(1-type1))-(1-type2) }

  root<- uniroot(f, nrange, theta=theta, sigma=sigma,
                type1=type1, type2=type2, l=l, u=u)$root
  n1<- ceiling(root)
  n2<- ceiling(root*r)
  n <- n1+n2
  return(list(n1=n1, n2=n2))
}
