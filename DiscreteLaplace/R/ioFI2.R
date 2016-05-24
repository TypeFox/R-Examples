ioFI2 <- function(x)
{
par <- estdlaplace2(x, "ML")
p <- par[1]
q <- par[2]
Dp2 <-
function(x,p,q)
{
sum(x<0)*(-1-log(p))/(p*log(p))^2+length(x)*(log(p*q)+1)/(p*log(p*q))^2-sum(x[x>=0])/p^2-sum(x>=0)/(1-p)^2
}
Dpq <-
function(x,p,q)
{
length(x)/(p*q*(log(p*q))^2)
}
Dq2 <-
function(x,p,q)
{
sum(x>=0)*(-1-log(q))/(q*log(q))^2+length(x)*(log(p*q)+1)/(q*log(p*q))^2+sum(x[x<0])/q^2+sum(x<0)*(1-2*q)/(q*(1-q))^2
}
solve(matrix(c(-Dp2(x,p,q),-Dpq(x,p,q),-Dpq(x,p,q),-Dq2(x,p,q)),2,2))
}
