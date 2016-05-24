"summary.grouped" <-
function(object, ...){
    if(!inherits(object, "grouped"))
        stop("Use only with 'grouped' objects.\n")
    dets <- object$details
    coefs <- object$coef
    p <- length(coefs)
    bname <- names(coefs)
    se <- sqrt(diag(vcov(object)))
    t.vals <- coefs[-p] / se[-p]
    p.vals <- 2 * ( 1 - pt(abs(t.vals), df = dets$n - p + 1))     
    coef.tab <- cbind(value = coefs[-p], st.err = se[-p], t.vals, p.vals)
    out <- list(object = object, coefficients = coef.tab, sigma = coefs[p], se.sigma = se[p],
            logLik = dets$logLik, AIC = AIC(object), BIC = AIC(object, k = log(dets$n)))
    class(out) <- "summ.grouped"
    out
}

