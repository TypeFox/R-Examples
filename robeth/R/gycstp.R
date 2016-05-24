"gycstp" <-
function(icase=.dFvGet()$ics,ialg=.dFvGet()$ilg,ni,a,e,tol=.dFvGet()$tlo,
maxit=.dFvGet()$mxt,t) {
if (missing(ni)) messagena("ni")
if (missing(a)) messagena("a")
if (missing(e)) messagena("e")
if (missing(t)) t <- single(1)
f.res <- .Fortran("gycstp",
icase=to.integer(icase),
ialg=to.integer(ialg),
ni=to.integer(ni),
a=to.single(a),
e=to.single(e),
tol=to.single(tol),
maxit=to.integer(maxit),
t=to.single(t))
list(t=f.res$t)
}
