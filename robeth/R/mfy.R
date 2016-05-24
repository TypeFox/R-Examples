"mfy" <-
function(a,y,m=nrow(a),iye=1,ize=1) {
if (missing(a)) messagena("a")
if (missing(y)) messagena("y")
n <- ncol(a)
mda <- nrow(a)
ny <- length(y)
nz <- ize*(m-1)+1
z <- single(nz)
f.res <- .Fortran("mfy",
a=to.single(a),
y=to.single(y),
z=to.single(z),
m=to.integer(m),
n=to.integer(n),
mda=to.integer(mda),
ny=to.integer(ny),
iye=to.integer(iye),
nz=to.integer(nz),
ize=to.integer(ize))
list(z=f.res$z)
}

