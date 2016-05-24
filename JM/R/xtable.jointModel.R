xtable.jointModel <-
function (x, caption = NULL, label = NULL, align = NULL, digits = NULL, 
    display = NULL, which = c("all", "Longitudinal", "Event"), varNames.Long = NULL, 
    varNames.Event = NULL, p.values = TRUE, digits.pval = 4, ...) {
    which <- match.arg(which)
    smobj <- summary(x)
    long <- smobj$'CoefTable-Long'
    event <- smobj$'CoefTable-Event'
    D <- smobj$D
    if (ncol(D) > 1 || nrow(D) == 1) {
        ind <- lower.tri(D, TRUE)
        D <- as.matrix(D[ind])
        rownames(D) <- paste("$D_{", 
            apply(which(ind, arr.ind = TRUE)[, 2:1, drop = FALSE], 1, 
                paste, sep = "", collapse = ""), "}$", sep = "")
        V <- vcov(x)
        ncolZ <- ncol(x$x$Z)
        L <- chol(x$coefficients$D)
        L <- t(L)[lower.tri(L, TRUE)]
        ii <- grep("B.", colnames(V), fixed = TRUE)
        V.re <- V[ii, ii]
        J <- jacobian2(L, ncolZ)
        se.D <- sqrt(diag(J %*% tcrossprod(V.re, J)))
    } else {
        V <- vcov(x)
        ii <- grep("B.", colnames(V), fixed = TRUE)
        V.re <- V[ii, ii]        
        se.D <- sqrt(diag(diag(c(D)) %*% tcrossprod(V.re, diag(c(D)))))
    }
    rnm.l <- c(rownames(long), "$\\log(\\sigma)$", NA, rownames(D))
    if (!is.null(varNames.Long) && is.character(varNames.Long) && 
        length(varNames.Long) == length(rnm.l))
        rnm.l <- varNames.Long
    se.logsigma <- sqrt(vcov(x)["Y.sigma", "Y.sigma"])
    pval.l <- c(long[, "p-value"], rep(NA, length(D) + 2))
    pval.l <- sprintf(paste("%.", digits.pval, "f", sep = ""), pval.l)
    if (any(ind <- pval.l == "0.0000"))
        pval.l[ind] <- "<0.0001"
    pval.l <- paste("$", pval.l, "$", sep = "")
    pval.l[pval.l == "$NA$"] <- as.character(NA)
    Dat.long <- data.frame(
        " " = rnm.l,
        Value = c(long[, "Value"], log(x$coefficients$sigma), NA, D),
        Std.Err = c(long[, "Std.Err"], se.logsigma, NA, se.D),
        "$p$-value" = pval.l,
        check.names = FALSE,
        row.names = seq_len(nrow(long) + 2 + length(D))
    )
    if (!p.values)
        Dat.long$"$p$-value" <- NULL
    rnm.e <- rownames(event)
    ii <- grep("log(xi.", rnm.e, fixed = TRUE)
    rnm.e[ii] <- paste("$\\log(\\xi_", seq_along(ii), ")$", sep = "")
    if (!is.null(varNames.Event) && is.character(varNames.Event) && 
        length(varNames.Event) == length(rnm.e))
        rnm.e <- varNames.Event
    pval.e <- event[, "p-value"]
    pval.e <- sprintf(paste("%.", digits.pval, "f", sep = ""), pval.e)
    if (any(ind <- pval.e == "0.0000"))
        pval.e[ind] <- "<0.0001"
    pval.e <- paste("$", pval.e, "$", sep = "")
    pval.e[pval.e == "$NA$"] <- as.character(NA)
    pval.e[ii] <- as.character(NA)
    Dat.surv <- data.frame(
        " " = rnm.e,
        Value = event[, "Value"],
        Std.Err = event[, "Std.Err"],
        "$p$-value" = pval.e,
        check.names = FALSE,
        row.names = seq_len(nrow(event))
    )
    if (!p.values)
        Dat.surv$"$p$-value" <- NULL
    nn <- nrow(Dat.long) - nrow(Dat.surv)
    Dat <- if (which == "all") {
        if (is.null(caption))
            caption <- paste("Parameter estimates, standard errors and $p$-values",
                "under the joint modeling analysis.",
                "$D_{ij}$ denote the $ij$-element of the covariance matrix",
                "for the random effects.")
        if (nn > 0) {
            dd <- Dat.surv[seq_len(nn), ]
            dd[] <- lapply(dd, function (x) {x[] <- NA; x})
            Dat.surv <- rbind(Dat.surv, dd)
        }
        if (nn < 0) {
            dd <- Dat.long[seq_len(abs(nn)), ]
            dd[] <- lapply(dd, function (x) {x[] <- NA; x})
            Dat.long <- rbind(Dat.long, dd)
        }
        lis <- if (p.values) {
            list(pos = list(-1), 
            command = paste("\\hline\n  & \\multicolumn{3}{c}{Event Process} &", 
                "& \\multicolumn{3}{c}{Longitudinal Process}\\\\\n"))
        } else {
            list(pos = list(-1), 
            command = paste("\\hline\n  & \\multicolumn{2}{c}{Event Process} &", 
                "& \\multicolumn{2}{c}{Longitudinal Process}\\\\\n"))
        }
        align <- if (p.values) "llrrrlrrr" else "llrrlrr"
        hline.after <- c(0, max(nrow(Dat.long), nrow(Dat.surv)))
        cbind(Dat.surv, Dat.long)
    } else if (which == "Longitudinal") {
        if (is.null(caption))
            caption <- paste("Parameter estimates, standard errors and $p$-values",
                "under the joint modeling analysis for the longitudinal linear mixed effects submodel.",
                "$D_{ij}$ denote the $ij$-element of the covariance matrix",
                "for the random effects.")
        lis <- NULL
        align <- if (p.values) "llrrr" else "llrr"
        hline.after <- c(-1, 0, nrow(Dat.long))
        Dat.long
    } else {
        if (is.null(caption))
            caption <- paste("Parameter estimates, standard errors and $p$-values",
                "under the joint modeling analysis for the event time submodel.")
        lis <- NULL
        align <- if (p.values) "llrrr" else "llrr"
        hline.after <- c(-1, 0, nrow(Dat.surv))
        Dat.surv
    }  
    if (!requireNamespace("xtable", quietly = TRUE)) 
        stop("'xtable' is required.\n")
    print(xtable(Dat, caption = caption, label = label, 
        align = align, digits = digits, display = display), 
        sanitize.text.function = function (x) x, include.rownames = FALSE,
        add.to.row = lis, hline.after = hline.after, ...)
}
