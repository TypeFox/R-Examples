d00.test <-
function (x) 
{

n <- length(x)
x.acf <- acf(x, plot=FALSE, lag.max=n)$acf[,,1]

l.acf <- length(x.acf)

mvec <- seq(from=1, to=l.acf, by=2)
svec <- (-1)^(2:(length(mvec)+1))
sel.acf <- x.acf[mvec+1]

T <- 4*sum(svec*sel.acf/mvec)/(sqrt(2)*pi)


polygamma.u <- 1.5 + (n-2)/2

varT <- (pi^2 - 2*psigamma(polygamma.u, deriv=1))/(n*pi^2)

Z <- T/sqrt(varT)

p.value <- 2*(1 - pnorm(abs(Z)))

#ll <- list(mvec=mvec, svec=svec, sel.acf=sel.acf, T=T, p.value=p.value,
#	statistic=Z)
ll <- list(p.value=p.value, statistic=Z, method="d00 test on acfs")

class(ll) <- "htest"

return(ll)
}
