"xsy" <-
function(x,y,s) {
if (missing(x)) messagena("x")
if (missing(y)) messagena("y")
if (missing(s)) messagena("s")
n <- length(x)
nn <- length(s)
result <- single(1)
f.res <- .Fortran("xsy",
x=to.single(x),
y=to.single(y),
s=to.single(s),
n=to.integer(n),
nn=to.integer(nn),
result=to.single(result))
list(result=f.res$result)
}

