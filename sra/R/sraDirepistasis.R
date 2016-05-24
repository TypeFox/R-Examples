sraDirepistasis <-
function (sradata, start = NULL, fixed = NULL, macroE = FALSE, 
    ...) 
{
    ss <- coef(sraCstvar(sradata = sradata))
    if (!is.null(start$mu0) && is.na(start$mu0)) {
        start$mu0 <- ss["mu0"]
    }
    if (!is.null(start$logvarA0) && is.na(start$logvarA0)) {
        start$logvarA0 <- ss["logvarA0"]
    }
    if (!is.null(start$logvarE0) && is.na(start$logvarE0)) {
        start$logvarE0 <- ss["logvarE0"]
    }
    default.start <- list(mu0 = NA, logvarA0 = NA, logvarE0 = NA, 
        logNe = NA)
    default.fixed <- list(logvarM = log(1e-20), logvarepsilon = log(1e-20))
    if (macroE) {
        default.start$logvarME <- NA
    }
    else {
        default.fixed$logvarME <- log(1e-20)
    }
    default.start[names(start)] <- start
    default.start[names(fixed)] <- NULL
    default.fixed[names(fixed)] <- fixed
    default.fixed[names(start)] <- NULL
    start <- default.start
    fixed <- default.fixed
    start[is.na(start)] <- sapply(names(start[is.na(start)]), 
        sraStartingvalues, sradata = sradata)
    nn <- names(start)
    start <- unlist(start)
    names(start) <- nn
    start = as.list(start)
    mlewrapper <- function(mu0, logvarA0, logvarE0, logNe, logvarM, 
        logepsilon, logminusepsilon, logvarepsilon, logvarME) {
        sraMinuslogL(sradata = sradata, FUNtimeseries = sraEpiTimeseries, 
            mu0 = mu0, logvarA0 = logvarA0, logvarE0 = logvarE0, 
            logNe = logNe, logepsilon = logepsilon, logminusepsilon = logminusepsilon, 
            logvarepsilon = logvarepsilon, logvarME = logvarME)
    }
    fit.pos <- try(mle(minuslogl = mlewrapper, start = c(start, 
        list(logepsilon = log(0.1))), fixed = c(fixed, list(logminusepsilon = NA)), 
        ...))
    fit.neg <- try(mle(minuslogl = mlewrapper, start = c(start, 
        list(logminusepsilon = log(0.1))), fixed = c(fixed, list(logepsilon = NA)), 
        ...))
    AIC.neg <- try(AIC(fit.neg), silent = TRUE)
    if (inherits(AIC.neg, "try-error")) {
        AIC.neg <- Inf
    }
    AIC.pos <- try(AIC(fit.pos), silent = TRUE)
    if (inherits(AIC.pos, "try-error")) {
        AIC.pos <- Inf
    }
    if (AIC.pos < AIC.neg) {
        fit <- fit.pos
        start <- c(start, list(logepsilon = coef(fit.pos)[["logepsilon"]]))
        fixed <- c(fixed, list(logminusepsilon = NA))
    }
    else {
        fit <- fit.neg
        start <- c(start, list(logminusepsilon = coef(fit.neg)[["logminusepsilon"]]))
        fixed <- c(fixed, list(logepsilon = NA))
    }
    return(sraMakeObject(sradata = sradata, model = fit, start = start, 
        fixed = fixed, FUNtimeseries = sraEpiTimeseries))
}
