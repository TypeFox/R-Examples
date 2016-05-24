"cerf" <-
function(x) {
if (missing(x)) messagena("x")
f <- single(1)
f.res <- .Fortran("cerf",
x=to.single(x),
f=to.single(f))
list(f=f.res$f)
}

