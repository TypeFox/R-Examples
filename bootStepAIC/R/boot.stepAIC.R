boot.stepAIC <-
function (object, data, B = 100, alpha = 0.05, direction = "backward", k = 2, 
                         verbose = FALSE, ...) {
    if (!class(object)[1] %in% c("lm", "aov", "glm", "negbin", "polr", "survreg", "coxph"))
        stop("\nboot.stepAIC() currently works for `lm', `aov', `glm', `negbin', `polr', `survreg' and `coxph' objects.\n")
    n <- nrow(data)
    index <- sample(n, n * B, replace = TRUE)
    dim(index) <- c(n, B)
    wrk.ind <- class(object)[1] %in% c("negbin", "survreg", "polr", "coxph")
    res <- vector(mode = "list", length = B)
    for (i in 1:B) {
        boot.data <- data[index[, i], ]
        try.newfit <- try({
            obj. <- if (inherits(object, "polr")) {
                update(object, data = boot.data, Hess = TRUE)
            } else { 
                update(object, data = boot.data)
            }
            if (wrk.ind) {
                Call <- obj.$call
                Call$data <- boot.data
                obj.$call <- Call
            }
            rr <- stepAIC(obj., direction = direction, trace = 0, k = k, ...)
        }, silent = TRUE)
        if (!inherits(try.newfit, "try-error")) {
            res[[i]] <- rr
            if (verbose)
                cat("\n", i, "replicate finished")
       }
    }
    # exclude fits that failed
    ind.fail <- sapply(res, is.null)
    if ((Bnew <- sum(!ind.fail)) < B)
        message("\nthe model fit failed in ", B - Bnew, " bootstrap samples")
    B <- Bnew
    res <- res[!ind.fail]
    # variable names
    vars <- lapply(res, function (x) if (length(tlabs <- attr(x$terms, "term.labels"))) tlabs else "Null")
    vars <- table(unlist(vars))
    nam.vars <- names(vars)
    Tab1 <- matrix(100 * as.vector(vars) / B, dimnames = list(nam.vars, "(%)"))
    ind.vars <- sapply(strsplit(nam.vars, ":"), length) < 2
    Tab1 <- Tab1[order(ind.vars, Tab1[, 1], decreasing = TRUE), , drop = FALSE]
    # coefficients sign
    coefs <- unlist(lapply(res, coef))
    coefs <- coefs[names(coefs) != "(Intercept)"]
    coefs <- lapply(split(coefs, names(coefs)), function (x) {
        nx <- length(x)
        c("+" = sum(x > 0) / nx, "-" = sum(x < 0) / nx)
    })
    nam.coefs <- names(coefs)
    Tab2 <- matrix(100 * unlist(coefs), ncol = 2, byrow = TRUE, dimnames = list(nam.coefs, c("+ (%)", "- (%)")))
    ind.coefs <- sapply(strsplit(nam.coefs, ":"), length) < 2
    Tab2 <- Tab2[order(ind.coefs, Tab2[, 1], decreasing = TRUE), , drop = FALSE]
    # statistical significance
    nams <- function (mat, name) {
        if (nrow(mat) == 1) {
            out <- mat[, name]
            names(out) <- rownames(mat)
            out
        } 
        else
            mat[, name]
    }
    pvals <- switch(class(object)[1],
        "lm" = unlist(lapply(res, function (x) nams(summary(x)$coef, "Pr(>|t|)"))),
        "aov" = unlist(lapply(res, function (x) {
            av <- anova(x)
            out <- av$"Pr(>F)"
            names(out) <- rownames(av)
            out
        })),
        "glm" = unlist(lapply(res, function (x) {
                sm <- summary(x) 
                if (sm$dispersion == 1) nams(sm$coef, "Pr(>|z|)") else nams(sm$coef, "Pr(>|t|)")
            })),
        "negbin" = unlist(lapply(res, function (x) nams(summary(x)$coef, "Pr(>|z|)"))),
        "polr" = 2 * pnorm(-abs(unlist(lapply(res, function (x) nams(summary(x)$coef, "t value"))))),
        "survreg" = unlist(lapply(res, function (x) nams(summary(x)$table, "p"))),
        "coxph" = unlist(lapply(res, function (x) if (!is.null(x$coef)) nams(summary(x)$coef, "Pr(>|z|)") else NA))
    )
    pvals <- pvals[!is.na(pvals)]
    pvals <- pvals[names(pvals) != "(Intercept)" & names(pvals) != "Log(scale)"]
    pvals <- lapply(split(pvals, names(pvals)), function (x) mean(x < alpha))
    nam.pvals <- names(pvals)
    Tab3 <- matrix(100 * unlist(pvals), ncol = 1, dimnames = list(nam.pvals, "(%)"))
    ind.pvals <- sapply(strsplit(nam.pvals, ":"), length) < 2
    Tab3 <- Tab3[order(ind.pvals, Tab3[, 1], decreasing = TRUE), , drop = FALSE]
    # Output
    out <- list("Covariates" = Tab1, "Sign" = Tab2, "Significance" = Tab3, "OrigModel" = object, 
                "OrigStepAIC" = stepAIC(object, direction = direction, trace = 0, k = k, ...), 
                direction = direction, k = k, "BootStepAIC" = res)
    class(out) <- "BootStep"
    out
}

