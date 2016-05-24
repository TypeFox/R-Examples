"print.summ.grouped" <-
function(x, ...){
    if(!inherits(x, "summ.grouped"))
        stop("Use only with 'summ.grouped' objects.\n")
    obj <- x$object
    dets <- obj$details
    cat("\nCall:\n", deparse(obj$call), "\n\n", sep = "")
    if(length(coefs <- x$coef)){
        cat("Model Summary:\n")
        model.sum <- data.frame("log.Lik" = round(x$logLik,3),
                                "AIC" = round(x$AIC,3),
                                "BIC" = round(x$BIC,3), row.names = "")
        print(model.sum)
        cat("\nCoefficients:\n")
        coefs <- data.frame("Esimate" = coefs[, 1], 
                            "Std.error" = coefs[, 2], 
                            "t.value" = coefs[, 3], 
                            "p.value" = format.pval(coefs[, 4], digits = 2, eps = 1e-03))
        if(nrow(coefs) == 1)
            row.names(coefs) <- "(Intercept)"
        print(coefs, digits = 3)
        cat("\nRandom-Effect:\n")
        distr <- dets$distr
        link.distr <- paste(dets$link, "-", distr, sep = "")
        distr <- if(distr == "t") paste(link.distr, "(df=", dets$df, ")", sep = "") else link.distr
        sigma <- data.frame("value" = x$sigma, 
                            "std.error" = x$se.sigma, 
                            "link-distribution" = distr, row.names = "sigma")
        print(sigma, digits = 3)
    } else cat("No coefficients\n")
    cat("\nOptimization:\n")
    cat("Convergence:", dets$conv, "\n")
    cat("max(|grad|):", format.pval(dets$max.sc, digits = 2, eps = 1e-06), "\n")
    cat(" Outer iter:", dets$k, "\n\n")
    invisible(x)
}

