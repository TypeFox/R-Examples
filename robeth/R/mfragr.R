"mfragr" <-
function(x,y,vp,nc,itype=.dFvGet()$ith,c=.dFvGet()$ccc,tol=.dFvGet()$tlo,gam=.dFvGet()$gma,
maxit=.dFvGet()$mxt,sigmac,sigmar) {
if (missing(y)) messagena("y")
if (missing(vp)) messagena("vp")
if (missing(nc)) messagena("nc")
if (missing(sigmac)) messagena("sigmac")
if (missing(sigmar)) messagena("sigmar")
n <- length(y)
np <- ncol(x)
mdx <- nrow(x)
ncov <- np*(np+1)/2
if (missing(x)) x <- matrix(single(1),mdx,np)
cpc <- single(nc)
cpr <- single(nc)
ipc <- integer(nc)
ipr <- integer(nc)
sc1 <- single(n)
sc2 <- single(n)
sc3 <- single(n)
sc4 <- single(n)
sc5 <- single(np)
sc6 <- single(np)
sc7 <- single(ncov)
ia <- integer(np)
ib <- integer(np)
ic <- integer(np)
f.res <- .Fortran("mfragr",
x=to.single(x),
y=to.single(y),
vp=to.single(vp),
n=to.integer(n),
np=to.integer(np),
mdx=to.integer(mdx),
ncov=to.integer(ncov),
nc=to.integer(nc),
itype=to.integer(itype),
c=to.single(c),
tol=to.single(tol),
gam=to.single(gam),
maxit=to.integer(maxit),
sigmac=to.single(sigmac),
sigmar=to.single(sigmar),
cpc=to.single(cpc),
cpr=to.single(cpr),
ipc=to.integer(ipc),
ipr=to.integer(ipr),
sc1=to.single(sc1),
sc2=to.single(sc2),
sc3=to.single(sc3),
sc4=to.single(sc4),
sc5=to.single(sc5),
sc6=to.single(sc6),
sc7=to.single(sc7),
ia=to.integer(ia),
ib=to.integer(ib),
ic=to.integer(ic))
list(x=f.res$x,cpc=f.res$cpc,cpr=f.res$cpr,ipc=f.res$ipc,ipr=f.res$ipr)
}
