"rybifr" <-
function(x,y,np,nthet=np+1,itype=2,icoll=0,isigma=1,ch=1.345,ck=
1.05*sqrt(nthet),bm=1.05*sqrt(nthet),tol=1.e-3,tau=1.e-6,maxitt=50,maxitw=80) {
if (missing(y)) messagena("y")
if (missing(np)) messagena("np")
n <- nrow(x)
ncov <- nthet*(nthet+1)/2
mdd <- nthet*(n+nthet+4)
mdr <- 2*(nthet+n)*(nthet+1)
mdi <- nthet
if (missing(x)) x <- matrix(single(1),n,nthet)
p <- ncol(x);if (nthet==p+1) x <- cbind(x,1) 
sigmaf <- single(1)
theta <- single(nthet)
rs <- single(n)
wgt <- single(n)
cov <- single(ncov)
dwrk <- double(mdd)
rwrk <- single(mdr)
iwrk <- integer(mdi)
f.res <- .Fortran("rybifr",
x=to.single(x),
y=to.single(y),
n=to.integer(n),
np=to.integer(np),
nthet=to.integer(nthet),
ncov=to.integer(ncov),
itype=to.integer(itype),
icoll=to.integer(icoll),
isigma=to.integer(isigma),
ch=to.single(ch),
ck=to.single(ck),
bm=to.single(bm),
tol=to.single(tol),
tau=to.single(tau),
maxitt=to.integer(maxitt),
maxitw=to.integer(maxitw),
sigmaf=to.single(sigmaf),
theta=to.single(theta),
rs=to.single(rs),
wgt=to.single(wgt),
cov=to.single(cov),
mdd=to.integer(mdd),
mdr=to.integer(mdr),
mdi=to.integer(mdi),
dwrk=as.double(dwrk),
rwrk=to.single(rwrk),
iwrk=to.integer(iwrk))
list(x=f.res$x,sigmaf=f.res$sigmaf,theta=f.res$theta,rs=f.res$rs,wgt=f.res$wgt,
cov=f.res$cov)
}

