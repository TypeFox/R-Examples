iFI2 <- function(p,q)
{
edp2 <- function(p,q)
{
log(q)/log(p)*(-(1-p)^2-p*log(p)*log(p*q))/((log(p*q))^2*p^2*(1-p)^2)
}
edpq <- function(p,q)
{
1/(p*q*(log(p*q))^2)
}
edq2 <- function(p,q)
{
log(p)/log(q)*(-(1-q)^2-q*log(q)*log(p*q))/((log(p*q))^2*q^2*(1-q)^2)
}
solve(matrix(c(-edp2(p,q),-edpq(p,q),-edpq(p,q),-edq2(p,q)),2,2))
}
