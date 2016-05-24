"kiedch" <-
function(wgt,c=.dFvGet()$ccc,itype=.dFvGet()$ite) {
if (missing(wgt)) messagena("wgt")
n <- length(wgt)
d <- single(n)
e <- single(n)
f.res <- .Fortran("kiedch",
wgt=to.single(wgt),
n=to.integer(n),
c=to.single(c),
itype=to.integer(itype),
d=to.single(d),
e=to.single(e))
list(d=f.res$d,e=f.res$e)
}
