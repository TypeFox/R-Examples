"tftaut" <-
function(rs1,rs2,wgt,exrho=rho,np,nq,sigma,itype=.dFvGet()$ite) {
if (missing(rs1)) messagena("rs1")
if (missing(rs2)) messagena("rs2")
if (missing(wgt)) messagena("wgt")
if (missing(np)) messagena("np")
if (missing(nq)) messagena("nq")
if (missing(sigma)) messagena("sigma")
n <- length(rs1)
sum1 <- single(1)
sum2 <- single(1)
ftau <- single(1)
f.res <- .Fortran("int52",
rs1=to.single(rs1),
rs2=to.single(rs2),
wgt=to.single(wgt),
as.integer(exrho()),
n=to.integer(n),
np=to.integer(np),
nq=to.integer(nq),
sigma=to.single(sigma),
itype=to.integer(itype),
sum1=to.single(sum1),
sum2=to.single(sum2),
ftau=to.single(ftau))
list(sum1=f.res$sum1,sum2=f.res$sum2,ftau=f.res$ftau)
}
