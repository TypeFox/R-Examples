"ruben" <-
function(xlmbda,delta,mult,x,xmode=1.0,maxit=50,eps=0.0001) {
if (missing(xlmbda)) messagena("xlmbda")
if (missing(delta)) messagena("delta")
if (missing(mult)) messagena("mult")
if (missing(x)) messagena("x")
n <- length(xlmbda)
dnsty <- single(1)
cumdf <- single(1)
ifault <- integer(1)
sg <- single(n)
st <- single(n)
sa <- single(maxit)
sb <- single(maxit)
f.res <- .Fortran("ruben",
xlmbda=to.single(xlmbda),
delta=to.single(delta),
mult=to.integer(mult),
n=to.integer(n),
x=to.single(x),
xmode=to.single(xmode),
maxit=to.integer(maxit),
eps=to.single(eps),
dnsty=to.single(dnsty),
cumdf=to.single(cumdf),
ifault=to.integer(ifault),
sg=to.single(sg),
st=to.single(st),
sa=to.single(sa),
sb=to.single(sb))
list(dnsty=f.res$dnsty,cumdf=f.res$cumdf,ifault=f.res$ifault)
}

