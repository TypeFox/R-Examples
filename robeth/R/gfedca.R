"gfedca" <-
function(vtheta,ci,wa,ni,oi=0,icase=.dFvGet()$ics) {
if (missing(vtheta)) messagena("vtheta")
if (missing(ci)) messagena("ci")
if (missing(wa)) messagena("wa")
if (missing(ni)) messagena("ni")
n <- length(vtheta)
d <- single(n)
e <- single(n)
if (length(oi)==1) oi <- rep(0,n)
f.res <- .Fortran("gfedca",
vtheta=to.single(vtheta),
ci=to.single(ci),
wa=to.single(wa),
ni=to.integer(ni),
oi=to.single(oi),
n=to.integer(n),
icase=to.integer(icase),
d=to.single(d),
e=to.single(e))
list(d=f.res$d,e=f.res$e)
}
