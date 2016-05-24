"Random" <-
function(iseed=.dFvGet()$isd) {
rn <- single(1)
f.res <- .Fortran("randow",
iseed=to.integer(iseed),
rn=to.single(rn))
z.res <- .Fortran("zdfvals",io=to.integer(0),dfv=single(66))
zdf <- z.res$dfv
zdf[48] <- to.single(f.res$iseed)
z.res <- .Fortran("zdfvals",io=to.integer(1),dfv=to.single(zdf))
f.res$rn
}
