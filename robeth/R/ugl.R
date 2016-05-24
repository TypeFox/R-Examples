"ugl" <-
function(upar,npar=4,s) {
if (missing(upar)) return (12)
if (missing(s)) messagena("s")
f.res <- .Fortran("int70",
upar=to.single(upar),
npar=to.integer(npar),
s=to.single(s),result=double(1))
f.res$result
}

