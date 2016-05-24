Edlaplace2 <-
function(p,q)
{
E1 <- log(q)/log(p*q)*p/(1-p)-log(p)/log(p*q)/(1-q)
E2 <- log(q)/log(p*q)*p*(1+p)/(1-p)^2+log(p)/log(p*q)*(1+q)/(1-q)^2
return(list(E1=E1,E2=E2))
}