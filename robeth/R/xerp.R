"xerp" <-
function(ip,xlcnst=-1.0,s) {
if (missing(ip)) messagena("ip")
if (missing(s)) messagena("s")
f <- single(1)
f.res <- .Fortran("xerp",
ip=to.integer(ip),
xlcnst=to.single(xlcnst),
s=to.single(s),
f=to.single(f))
list(xlcnst=f.res$xlcnst,f=f.res$f)
}

