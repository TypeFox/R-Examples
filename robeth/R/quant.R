"quant" <-
function(p) {
if (missing(p)) messagena("p")
x <- single(1)
f.res <- .Fortran("nquant",
p=to.single(p),
x=to.single(x))
list(x=f.res$x)
}

