"lytau2" <-
function(z,expsi=psi,expsp=psp,exchi=chi,exrho=rho,m,n,tol=.dFvGet()$tlo,
gam=.dFvGet()$gma,isigma=.dFvGet()$isg,maxit=.dFvGet()$mxt,nitmon=.dFvGet()$ntm) {
if (missing(z)) messagena("z")
if (missing(m)) messagena("m")
if (missing(n)) messagena("n")
mpn <- length(z)
thetal <- single(1)
deltal <- single(1)
thetas <- single(1)
sigmaf <- single(1)
ftau <- single(1)
p <- single(1)
rs1 <- single(mpn)
rs2 <- single(mpn)
cov <- single(3)
work1 <- matrix(single(1),mpn,6)
work2 <- single(8)
iwork <- integer(2)
f.res <- .Fortran("int36",
z=to.single(z),
as.integer(expsi()),
as.integer(expsp()),
as.integer(exchi()),
as.integer(exrho()),
m=to.integer(m),
n=to.integer(n),
mpn=to.integer(mpn),
tol=to.single(tol),
gam=to.single(gam),
isigma=to.integer(isigma),
maxit=to.integer(maxit),
nitmon=to.integer(nitmon),
thetal=to.single(thetal),
deltal=to.single(deltal),
thetas=to.single(thetas),
sigmaf=to.single(sigmaf),
ftau=to.single(ftau),
p=to.single(p),
rs1=to.single(rs1),
rs2=to.single(rs2),
cov=to.single(cov),
work1=to.single(work1),
work2=to.single(work2),
iwork=to.integer(iwork))
list(thetal=f.res$thetal,deltal=f.res$deltal,thetas=f.res$thetas,sigmaf=
f.res$sigmaf,ftau=f.res$ftau,p=f.res$p,rs1=f.res$rs1,rs2=f.res$rs2,cov=f.res$cov)
}
