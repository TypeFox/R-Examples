"dpoiss" <-
function(y,ci,vtheta,wa,f0,oi=0,kap=1e-6) {
if (missing(y)) messagena("y")
if (missing(ci)) messagena("ci")
if (missing(vtheta)) messagena("vtheta")
if (missing(wa)) messagena("wa")
if (missing(f0)) messagena("f0")
n <- length(y)
d <- single(n)
if(length(oi)==1) oi <- rep(0,n)
f.res <- .Fortran("dpoiss",
y=to.single(y),
ci=to.single(ci),
vtheta=to.single(vtheta),
wa=to.single(wa),
f0=to.single(f0),
oi=to.single(oi),
n=to.integer(n),
kap=to.single(kap),
d=to.single(d))
list(d=f.res$d)
}

