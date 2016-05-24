summary.LORgee <-
function (object , ...) 
{
    robust.se <- sqrt(diag(object$robust.variance))
    robust.z <- coef(object)/robust.se
    pvalue <- 2 * (1 - pnorm(abs(robust.z)))
    TAB <- cbind(Estimate = coef(object), Robust_S.E. = robust.se, 
        Robust_z = robust.z, p.value = pvalue)
    TAB <- round(TAB, 5)
    colnames(TAB) <- c("Estimate", "san.se", "san.z", "Pr(>|san.z|)")
    res <- list(title = object$title, version = object$version, coefficients = TAB, 
        link = object$link, local.odds.ratios = object$local.odds.ratios, call = object$call, 
        residuals = object$residuals, niter = object$convergence$niter, 
        pvalue = object$pvalue)
    class(res) <- "summary.LORgee"
    res
}

