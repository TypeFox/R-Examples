### Since we had an undetected bug in  rnorMix()...

## These are defined as norMix() calls in  ../R/zMarrWand-dens.R
library("nor1mix")
ppos <- which("package:nor1mix" == search())
nms <- ls(pat="^MW.nm", pos = ppos)
nms <- nms[order(as.numeric(substring(nms,6)))] # warning <== "MW.nm2.old"

set.seed(123)
for(n in nms) {
    o <- get(n, pos = ppos)
    cat("\n",n,":\n"); print(o)
    cat("4 random X from", n,":")
    print(rnorMix(4, o))

    ## Testing of sort.norMix():
    if(is.unsorted(o[,"mu"]))
        o <- sort(o)
    ae <- all.equal(o, sort(o), tolerance = 4e-16) # non-stable sorting, MW.
    if(!isTRUE(ae)) {
        cat("sort(o) differs from o:\n"); print(ae)
    }
}
