sraAutoregEvolv <-
function (sradata, active = c(FALSE, TRUE, FALSE, FALSE), start = NULL, 
    fixed = NULL, negative.k = FALSE, rand = 0, rep = 1, ...) 
{
    if (rep > 1) {
        ans <- NULL
        if (rand == 0) {
            rand <- 0.01
        }
        for (r in 1:rep) {
            fit <- try(sraAutoregEvolv(sradata = sradata, active = active, 
                start = start, fixed = fixed, negative.k = negative.k, 
                rand = rand, rep = 1, ...))
            if (!inherits(fit, "try-error") && (is.null(ans) || 
                (AIC(fit) < AIC(ans)))) {
                ans <- fit
            }
        }
        return(ans)
    }
    default.start <- list(mu0 = NA, logIA0 = NA, logIE0 = NA)
    if (negative.k) {
        if (active[1]) {
            default.start$relativekA0 = NA
            default.start$relativekE0 = NA
        }
        if (active[2]) {
            default.start$kA1 = NA
            default.start$kE1 = NA
        }
        if (active[3]) {
            default.start$kA2 = NA
            default.start$kE2 = NA
        }
        if (active[4]) {
            default.start$kA3 = NA
            default.start$kE3 = NA
        }
    }
    else {
        if (active[1]) {
            default.start$relativekA0 = NA
            default.start$relativekE0 = NA
        }
        if (active[2]) {
            default.start$logkA1 = NA
            default.start$logkE1 = NA
        }
        if (active[3]) {
            default.start$logkA2 = NA
            default.start$logkE2 = NA
        }
        if (active[4]) {
            default.start$logkA3 = NA
            default.start$logkE3 = NA
        }
    }
    default.start[names(start)] <- start
    start <- default.start
    start[names(fixed)] <- NULL
    default.fixed <- list(mu0 = 0, logIA0 = 0, logIE0 = 0, relativekA0 = 0, 
        relativekE0 = 0, kA1 = 1, kE1 = 1, kA2 = 0, kE2 = 0, 
        kA3 = 0, kE3 = 0)
    if (!negative.k) {
        default.fixed <- list(mu0 = 0, logIA0 = 0, logIE0 = 0, 
            relativekA0 = 0, relativekE0 = 0, logkA1 = 0, logkE1 = 0, 
            logkA2 = -20, logkE2 = -20, logkA3 = -20, logkE3 = -20)
    }
    default.fixed[names(fixed)] <- fixed
    default.fixed[names(start)] <- NULL
    fixed <- default.fixed
    start[is.na(start)] <- sapply(names(start[is.na(start)]), 
        sraStartingvalues, sradata = sradata)
    mlewrapper.neg <- function(mu0, logIA0, logIE0, relativekA0, 
        kA1, kA2, kA3, relativekE0, kE1, kE2, kE3) {
        sraMinuslogL(sradata = sradata, FUNtimeseries = sraAutoregEvolvTimeseries, 
            mu0 = mu0, logIA0 = logIA0, logIE0 = logIE0, relativekA0 = relativekA0, 
            kA1 = kA1, kA2 = kA2, kA3 = kA3, relativekE0 = relativekE0, 
            kE1 = kE1, kE2 = kE2, kE3 = kE3)
    }
    mlewrapper.noneg <- function(mu0, logIA0, logIE0, relativekA0, 
        logkA1, logkA2, logkA3, relativekE0, logkE1, logkE2, 
        logkE3) {
        sraMinuslogL(sradata = sradata, FUNtimeseries = sraAutoregEvolvTimeseries, 
            mu0 = mu0, logIA0 = logIA0, logIE0 = logIE0, relativekA0 = relativekA0, 
            kA1 = exp(logkA1), kA2 = exp(logkA2), kA3 = exp(logkA3), 
            relativekE0 = relativekE0, kE1 = exp(logkE1), kE2 = exp(logkE2), 
            kE3 = exp(logkE3))
    }
    fit <- NULL
    if (negative.k) {
        fit <- mle(minuslogl = mlewrapper.neg, start = start, 
            fixed = fixed, ...)
    }
    else {
        fit <- mle(minuslogl = mlewrapper.noneg, start = start, 
            fixed = fixed, ...)
    }
    return(sraMakeObject(sradata = sradata, model = fit, start = start, 
        fixed = fixed, FUNtimeseries = sraAutoregEvolvTimeseries))
}
