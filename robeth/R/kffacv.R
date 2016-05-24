"kffacv" <-
function(rs,expsi=psi,expsp=psp,np,sigma) {
if (missing(rs)) messagena("rs")
if (missing(np)) messagena("np")
if (missing(sigma)) messagena("sigma")
n <- length(rs)
fh <- single(1)
f.res <- .Fortran("int25",
rs=to.single(rs),
as.integer(expsi()),
as.integer(expsp()),
n=to.integer(n),
np=to.integer(np),
sigma=to.single(sigma),
fh=to.single(fh))
list(fh=f.res$fh)
}

