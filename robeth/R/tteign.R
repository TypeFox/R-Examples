"tteign" <-
function(covtau,nq,mdc=np-nq) {
if (missing(covtau)) messagena("covtau")
if (missing(nq)) messagena("nq")
np <- ncol(covtau)
xlmbda <- single(np)
iv <- integer(np)
sv <- single(np)
f.res <- .Fortran("tteign",
covtau=to.single(covtau),
np=to.integer(np),
nq=to.integer(nq),
mdc=to.integer(mdc),
xlmbda=to.single(xlmbda),
iv=to.integer(iv),
sv=to.single(sv))
list(xlmbda=f.res$xlmbda)
}

