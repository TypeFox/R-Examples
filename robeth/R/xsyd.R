"xsyd" <-
function(x,y,s) {
if (missing(x)) messagena("x")
if (missing(y)) messagena("y")
if (missing(s)) messagena("s")
n <- length(x)
nn <- length(s)
result <- double(1)
f.res <- .Fortran("xsyd",
x=as.double(x),
y=as.double(y),
s=as.double(s),
n=to.integer(n),
nn=to.integer(nn),
result=as.double(result))
list(result=f.res$result)
}

