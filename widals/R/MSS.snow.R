MSS.snow <-
function (FUN.source, current.best, p.ndx.ls, f.d, sds.mx, k.glob, k.loc.coef, X = NULL)
{
    
    xenvr <- as.environment(1)
    envmh <- environment(NULL)
    GP <- GP
    if(is.function(FUN.source)) {
        FUN.source()
    } else {
        if (!is.null(FUN.source)) {
            source(FUN.source)
        }
    }
    if (is.na(current.best)) {
        GP.mx <- matrix(GP, 1, length(GP))
        if (!is.null(FUN.GP)) {
            GP.mx <- FUN.GP(GP.mx)
        }
        current.best <- FUN.MH(1, GP.mx = GP.mx, X = X)
        cat(current.best, "\n")
    }
    if (!is.null(k.glob)) {
        for (k.times in 1:k.glob) {
            cat(k.times, "\n")
            for (p.ndx in p.ndx.ls) {
                n.mh <- as.integer(k.loc.coef * 2^length(p.ndx))
                GP.mx <- matrix(GP, n.mh, length(GP), byrow = TRUE)
                for (ip in p.ndx) {
                    GP.mx[, ip] <- f.d[[ip]](n.mh, GP[ip], sds.mx[k.times,
                    ip])
                }
                if (!is.null(FUN.GP)) {
                    GP.mx <- FUN.GP(GP.mx)
                }
                if (run.parallel) {
                    sfOut <- sfClusterApplyLB(1:n.mh, FUN.MH, GP.mx = GP.mx,
                    X = X)
                }
                else {
                    sfOut <- list()
                    for (jj in 1:n.mh) {
                        sfOut[[jj]] <- FUN.MH(jj, GP.mx = GP.mx,
                        X = X)
                    }
                }
                errs <- unlist(sfOut)
                errs[is.na(errs)] <- Inf
                errs[is.nan(errs)] <- Inf
                best.ndx <- which(errs == min(errs))[1]
                if (errs[best.ndx] < current.best) {
                    current.best <- errs[best.ndx]
                    GP <- GP.mx[best.ndx, , drop = TRUE]
                    assign("current.best", current.best, envir = envmh)
                    assign("current.best.GP", GP, envir = envmh)
                    X <- FUN.I(envmh = envmh, X = X)
                }
            }
        }
    }
    if (!is.null(FUN.EXIT)) {
        FUN.EXIT(envmh = envmh, X = X)
    }
    assign("GP", GP, pos = xenvr)
}
