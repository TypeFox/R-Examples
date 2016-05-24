"exchd" <-
function(s,n,h,k) {
if (missing(n)) messagena("n")
if (missing(h)) messagena("h")
if (missing(k)) messagena("k")
nn <- length(s)
if (missing(s)) s <- double(nn)
f.res <- .Fortran("exchd",
s=as.double(s),
n=to.integer(n),
nn=to.integer(nn),
h=to.integer(h),
k=to.integer(k))
list(s=f.res$s)
}

