"libeth" <-
function(d=.dFvGet()$ddd) {
bta <- single(1)
f.res <- .Fortran("libeth",
d=to.single(d),
bta=to.single(bta))
list(bta=f.res$bta)
}
