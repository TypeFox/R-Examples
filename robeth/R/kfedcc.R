"kfedcc" <-
function(wgt,rs,expsi=psi,expsp=psp,sigma,itype=.dFvGet()$ite) {
if (missing(wgt)) messagena("wgt")
if (missing(rs)) messagena("rs")
if (missing(sigma)) messagena("sigma")
n <- length(wgt)
d <- single(n)
e <- single(n)
f.res <- .Fortran("int27",
wgt=to.single(wgt),
rs=to.single(rs),
as.integer(expsi()),
as.integer(expsp()),
n=to.integer(n),
sigma=to.single(sigma),
itype=to.integer(itype),
d=to.single(d),
e=to.single(e))
list(d=f.res$d,e=f.res$e)
}
