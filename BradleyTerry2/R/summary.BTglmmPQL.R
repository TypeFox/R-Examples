summary.BTglmmPQL <- function(object, dispersion = NULL, correlation = FALSE,
                            symbolic.cor = FALSE, ...) {
    if (identical(object$sigma, 0)){
        cat("PQL algorithm converged to fixed effects model\n")
        return(NextMethod("summary"))
    }
    aliased <- is.na(coef(object))
    coefs <- coef(object)[!aliased]
    cov.scaled <- cov.unscaled <- object$varFix # when dispersion != 1?
    dn <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    if (object$rank > 0) {
        sterr <- sqrt(diag(cov.scaled))
        tvalue <- coefs/sterr
        pvalue <- 2 * pnorm(-abs(tvalue))
        fixef.table <- cbind(coefs, sterr, tvalue, pvalue)
        dimnames(fixef.table) <- list(names(coefs), dn)
    }
    else {
        fixef.table <- matrix(, 0, 4)
        dimnames(fixef.table) <- list(NULL, dn)
    }
    sterr <- sqrt(object$varSigma)
    tvalue <- object$sigma/sterr
    pvalue <- 2 * pnorm(-abs(tvalue))
    ranef.table <- cbind(object$sigma, sterr, tvalue, pvalue)
    dimnames(ranef.table) <- list("Std. Dev.", dn)
    ans <- c(object[c("call", "family", "iter", "rank", "na.action")],
             list(fixef = fixef.table, ranef = ranef.table,
                  aliased = aliased, dispersion = 1,
                  cov.unscaled = cov.unscaled))
    if (correlation & object$rank > 0) {
        dd <- sqrt(diag(cov.unscaled))
        ans$correlation <- cov.unscaled/outer(dd, dd)
        ans$symbolic.cor <- symbolic.cor
    }
    class(ans) <- "summary.BTglmmPQL"
    ans
}
