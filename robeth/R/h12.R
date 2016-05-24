"h12" <-
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
if (missing(u)) u <- matrix(single(1),iue,m)
if (missing(up)) up <- single(1)
if (missing(c)) c <- single(mdc)
f.res <- .Fortran("h12",
mode=to.integer(mode),
lpivot=to.integer(lpivot),
l1=to.integer(l1),
m=to.integer(m),
u=to.single(u),
iue=to.integer(iue),
up=to.single(up),
c=to.single(c),
ice=to.integer(ice),
icv=to.integer(icv),
ncv=to.integer(ncv),
mdc=to.integer(mdc))
list(u=f.res$u,up=f.res$up,c=f.res$c)
}

