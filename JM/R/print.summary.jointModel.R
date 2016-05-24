print.summary.jointModel <-
function (x, digits = max(4, getOption("digits") - 4), 
        printKnots = FALSE, ...) {
    if (!inherits(x, "summary.jointModel"))
        stop("Use only with 'summary.jointModel' objects.\n")
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", 
        collapse = "\n"), "\n\n", sep = "")
    cat("Data Descriptives:\n")
    pcEv <- round(100 * sum(x$d) / x$n, 1)
    cat("Longitudinal Process\t\tEvent Process")
    cat("\nNumber of Observations: ", x$N, "\tNumber of Events: ", 
        sum(x$d), " (", pcEv, "%)", sep = "")
    cat("\nNumber of Groups:", length(unique(x$id)))
    cat("\n\nJoint Model Summary:")
    cat("\nLongitudinal Process: Linear mixed-effects model")
    cat("\nEvent Process: ")
    if (x$method == "Cox-PH-GH") {
        cat("Relative risk model with unspecified baseline risk function\n\n")
    } else if (x$method == "weibull-AFT-GH") {
        cat("Weibull accelerated failure time model\n")
    } else if (x$method == "weibull-PH-GH") {
        cat("Weibull relative risk model\n")
    } else if (x$method == "piecewise-PH-GH") {
        if (printKnots)
            cat("Relative risk model with piecewise-constant\n\tbaseline risk function (knots at: ", 
                paste(round(x$control$knots, 1), collapse = ", "), ")\n", sep = "")
        else 
            cat("Relative risk model with piecewise-constant\n\t\tbaseline risk function\n")
    } else if (x$method == "spline-PH-GH") {
        xx <- if (length(x$control$knots) == 1) {
            kk <- round(unique(x$control$knots[[1]]), 1)
            paste(kk[-c(1, length(kk))], collapse = ", ")
        } else {
            paste(names(x$control$knots), sapply(x$control$knots, function (k) {
                kk <- round(unique(k), 1)
                paste(kk[-c(1, length(kk))], collapse = ", ")
            }), sep = ": ", collapse = "\n\t\t")
        }
        if (printKnots)
            cat("Relative risk model with spline-approximated baseline risk function (knots at: ", xx, ")\n", sep = "")
        else {
            if (x$CompRisk) {
                cat("Competing risks relative risk model with spline-approximated\n\t\tbaseline risk function\n")
            } else if (x$nstrata > 1) {
                cat("Stratified relative risk model with spline-approximated\n\t\tbaseline risk function\n")
            } else {
                cat("Relative risk model with spline-approximated\n\t\tbaseline risk function\n")
            }   
        }
    } else {
        cat("log cumulative baseline hazard with B-splines (internal knots at: ", 
            paste(round(exp(x$knots[-c(1, length(x$knots))]), 2), collapse = ", "), ")\n", sep = "")
    }
    cat("Parameterization:", switch(x$parameterization, "value" = "Time-dependent", 
        "slope" = "Time-dependent slope", "both" = "Time-dependent + time-dependent slope"), "\n\n")
    model.sum <- data.frame(log.Lik = x$logLik, AIC = x$AIC, BIC = x$BIC, row.names = "")
    print(model.sum)
    cat("\nVariance Components:\n")
    D <- x$D
    ncz <- nrow(D)
    diag.D <- ncz != ncol(D)
    sds <- if (diag.D) sqrt(D) else sqrt(diag(D))
    if (ncz > 1) {
        if (diag.D) {
            dat <- as.data.frame(round(rbind(sds, "Residual" = x$sigma), digits))
            names(dat) <- "StdDev"
        } else {
            corrs <- cov2cor(D)
            corrs[upper.tri(corrs, TRUE)] <- 0
            mat <- round(cbind(sds, corrs[, -ncz]), digits)
            mat <- rbind(mat, c(x$sigma, rep(0, ncz - 1)))
            mat <- apply(mat, 2, sprintf, fmt = "% .4f")
            mat[mat == mat[1, 2]] <- ""
            mat[1, -1] <- abbreviate(colnames(D)[-ncz], 6)
            colnames(mat) <- c(colnames(mat)[1], rep("", ncz - 1))
            dat <- data.frame(mat, check.rows = FALSE, check.names = FALSE)
            names(dat) <- c("StdDev", "Corr", if (ncz > 2) rep(" ", ncz - 2) else NULL)
            row.names(dat) <- c(dimnames(D)[[1]], "Residual")
        }
    } else {
        dat <- data.frame("StdDev" = c(sds, x$sigma), row.names = c(rownames(D), "Residual"), 
            check.rows = FALSE, check.names = FALSE)
    }
    print(dat)
    cat("\nCoefficients:")
    cat("\nLongitudinal Process\n")
    out <- as.data.frame(round(x$"CoefTable-Long", digits))
    ind <- out$"p-value" == 0
    out$"p-value" <- sprintf(paste("%.", digits, "f", sep = ""), out$"p-value")
    out$"p-value"[ind] <- paste("<0.", paste(rep("0", digits - 1), collapse = ""), "1", sep = "")
    print(out)
    cat("\nEvent Process\n")
    out <- as.data.frame(round(x$"CoefTable-Event", digits))
    ind <- out$"p-value" == 0
    out$"p-value" <- sprintf(paste("%.", digits, "f", sep = ""), out$"p-value")
    out$"p-value"[ind] <- paste("<0.", paste(rep("0", digits - 1), collapse = ""), "1", sep = "")
    if (x$method == "piecewise-PH-GH") {
        Q <- length(x$control$knots) + 1
        ind.xi <- row.names(out) %in% paste("log(xi.", 1:Q, ")", sep = "")
        out$"p-value"[ind.xi] <- " "
    }
    print(out)
    if(x$method == "weibull-PH-GH" || x$method == "weibull-AFT-GH")
        cat("\nScale:", round(exp(x$"CoefTable-Event"[nrow(x$"CoefTable-Event"), 1]), digits), "\n")
    cat("\nIntegration:\n")
    GH <- x$method %in% c("ch-GH", "Cox-PH-GH", "weibull-PH-GH", "weibull-AFT-GH", "piecewise-PH-GH", "spline-PH-GH")
    cat("method:", if (GH && x$control$typeGH == "simple") "Gauss-Hermite" 
        else if (GH && x$control$typeGH != "simple") "(pseudo) adaptive Gauss-Hermite" else "Laplace")
    if (GH) cat("\nquadrature points:", x$control$GHk, "\n") else cat("\n")
    cat("\nOptimization:\n")
    cat("Convergence:", as.numeric(x$conv), "\n")
    cat("\n")
    invisible(x)    
}
