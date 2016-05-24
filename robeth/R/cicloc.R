"cicloc" <-
function(eps=.dFvGet()$esp,tol=.dFvGet()$tlo) {
c <- single(1)
f.res <- .Fortran("cicloc",
eps=to.single(eps),
tol=to.single(tol),
c=to.single(c))
list(c=f.res$c)
}
