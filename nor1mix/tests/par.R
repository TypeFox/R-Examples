library(nor1mix)

## Check  nM2par(), par2norMix() and llnorMix() :

nms <- paste("MW.nm", 1:16, sep="")
for(n in nms) {
    cat(n,":")
    obj <- get(n, envir = as.environment("package:nor1mix"))
    pp <- nM2par(obj)

    stopifnot(all.equal(pp, nM2par(par2norMix(pp)), tol= 1e-15),
              all.equal(obj, par2norMix(nM2par(obj)),
                        check.attributes=FALSE, tol=1e-15))
    xx <- rnorMix(1000, obj)
    stopifnot(all.equal(llnorMix(nM2par(obj), xx),
                        sum(dnorMix(xx, obj, log = TRUE)),
                    tol = 1e-15))
    cat(" [ok]\n")
}
