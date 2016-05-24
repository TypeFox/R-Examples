"wfshat" <-
function(xt,n=nrow(xt),sh) {
if (missing(xt)) messagena("xt")
if (missing(sh)) messagena("sh")
np <- ncol(xt)
mdx <- nrow(xt)
wgt <- single(n)
sc <- single(n)
f.res <- .Fortran("wfshat",
xt=to.single(xt),
n=to.integer(n),
np=to.integer(np),
mdx=to.integer(mdx),
wgt=to.single(wgt),
sh=to.single(sh),
sc=to.single(sc))
list(wgt=f.res$wgt)
}

