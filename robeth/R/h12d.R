"h12d" <-
function(mode,lpivot,l1,u,up,c,ice,icv,ncv) {
if (missing(mode)) messagena("mode")
if (missing(lpivot)) messagena("lpivot")
if (missing(l1)) messagena("l1")
if (missing(ice)) messagena("ice")
if (missing(icv)) messagena("icv")
if (missing(ncv)) messagena("ncv")
m <- ncol(u)
iue <- nrow(u)
mdc <- length(c)
if (missing(u)) u <- matrix(double(1),iue,m)
if (missing(up)) up <- double(1)
if (missing(c)) c <- double(mdc)
f.res <- .Fortran("h12d",
mode=to.integer(mode),
lpivot=to.integer(lpivot),
l1=to.integer(l1),
m=to.integer(m),
u=as.double(u),
iue=to.integer(iue),
up=as.double(up),
c=as.double(c),
ice=to.integer(ice),
icv=to.integer(icv),
ncv=to.integer(ncv),
mdc=to.integer(mdc))
list(u=f.res$u,up=f.res$up,c=f.res$c)
}

