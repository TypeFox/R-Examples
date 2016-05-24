"LSWsim"<- 
function(spec){
        #
        #
        # First check that all spectral elements are non-zero
        #
        if (any(spec$D < 0))
                stop("All spectral elements must be non-negative")
        #
        #
        # Now multiply by random element and factor of 2 (to undo AvBasis
        # averaging)
        #
        nlev <- nlevelsWT(spec)
        len <- 2^nlev
        for(i in (nlev-1):0)    {
                v <- accessD(spec, level=i)
                v <- sqrt(v)*2^(nlev-i)*rnorm(len)
                spec <- putD(spec, level=i, v=v)
                }
        AvBasis(convert(spec))
}
"cns"<- 
function(n, filter.number=1, family="DaubExPhase"){
        if (is.na(IsPowerOfTwo(n)))
                stop("n must be a power of two")
        z <- rep(0, n)
        zwdS <- wd(z, filter.number=filter.number, family=family, type="station")
        zwdS
        }
"checkmyews" <- function(spec, nsim=10){
        ans <- cns(2^nlevelsWT(spec))
        for(i in 1:nsim)        {
                cat(".")
                LSWproc <- LSWsim(spec)
                ews <- ewspec(LSWproc, filter.number=1, family="DaubExPhase",
                        WPsmooth=FALSE)
                ans$D <- ans$D + ews$S$D
                ans$C <- ans$C + ews$S$C
        }
        ans$D <- ans$D/nsim
        ans$C <- ans$C/nsim
        ans
}
