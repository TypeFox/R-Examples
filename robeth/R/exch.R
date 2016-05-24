"exch" <-
function(s,n,h,k) {
if (missing(n)) messagena("n")
if (missing(h)) messagena("h")
if (missing(k)) messagena("k")
nn <- length(s)
if (missing(s)) s <- single(nn)
f.res <- .Fortran("exch",
s=to.single(s),
n=to.integer(n),
nn=to.integer(nn),
h=to.integer(h),
k=to.integer(k))
list(s=f.res$s)
}

