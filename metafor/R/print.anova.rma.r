print.anova.rma <-
function (x, digits, ...) 
{
    if (class(x) != "anova.rma") 
        stop("Argument 'x' must be an object of class \"anova.rma\".")
    if (missing(digits)) 
        digits <- x$digits
    if (x$test == "Wald.b") {
        cat("\n")
        cat("Test of Moderators (coefficient(s) ", paste(x$btt, 
            collapse = ","), "): \n", sep = "")
        if (x$knha) {
            cat("F(df1 = ", x$m, ", df2 = ", x$dfs, ") = ", formatC(x$QM, 
                digits = digits, format = "f"), ", p-val ", .pval(x$QMp, 
                digits = digits, showeq = TRUE, sep = " "), "\n\n", 
                sep = "")
        }
        else {
            cat("QM(df = ", x$m, ") = ", formatC(x$QM, digits = digits, 
                format = "f"), ", p-val ", .pval(x$QMp, digits = digits, 
                showeq = TRUE, sep = " "), "\n\n", sep = "")
        }
    }
    if (x$test == "Wald.L") {
        cat("\n")
        if (x$m == 1) {
            cat("Hypothesis:")
        }
        else {
            cat("Hypotheses:")
        }
        print(x$hyp)
        cat("\nResults:\n")
        res.table <- cbind(estimate = x$Lb, se = x$se, zval = x$zval, 
            pval = x$pval)
        colnames(res.table) <- c("estimate", "se", "zval", "pval")
        if (x$knha) 
            colnames(res.table)[3] <- "tval"
        rownames(res.table) <- paste0(seq_len(x$m), ":")
        res.table <- formatC(res.table, digits = digits, format = "f")
        res.table[, 4] <- .pval(x$pval, digits = digits)
        print(res.table, quote = FALSE, right = TRUE)
        cat("\n")
        if (!is.null(x$QM)) {
            if (x$m == 1) {
                cat("Test of Hypothesis:\n")
            }
            else {
                cat("Omnibus Test of Hypotheses:\n")
            }
            if (x$knha) {
                cat("F(df1 = ", x$m, ", df2 = ", x$dfs, ") = ", 
                  formatC(x$QM, digits = digits, format = "f"), 
                  ", p-val ", .pval(x$QMp, digits = digits, showeq = TRUE, 
                    sep = " "), "\n\n", sep = "")
            }
            else {
                cat("QM(df = ", x$m, ") = ", formatC(x$QM, digits = digits, 
                  format = "f"), ", p-val ", .pval(x$QMp, digits = digits, 
                  showeq = TRUE, sep = " "), "\n\n", sep = "")
            }
        }
    }
    if (x$test == "LRT") {
        res.table <- rbind(c(x$p.f, x$fit.stats.f[3], x$fit.stats.f[4], 
            x$fit.stats.f[5], x$fit.stats.f[1], NA, NA, x$QE.f, 
            x$tau2.f, NA), c(x$p.r, x$fit.stats.r[3], x$fit.stats.r[4], 
            x$fit.stats.r[5], x$fit.stats.r[1], x$LRT, x$pval, 
            x$QE.r, x$tau2.r, NA))
        res.table[, seq.int(from = 2, to = 10)] <- formatC(res.table[, 
            seq.int(from = 2, to = 10)], digits = digits, format = "f")
        colnames(res.table) <- c("df", "AIC", "BIC", "AICc", 
            "logLik", "LRT", "pval", "QE", "tau^2", "R^2")
        rownames(res.table) <- c("Full", "Reduced")
        res.table[2, 7] <- .pval(x$pval, digits = digits)
        res.table[1, c(6, 7)] <- ""
        res.table[1, 10] <- ""
        res.table[2, 10] <- paste0(x$R2, "%")
        if (x$method == "FE" || is.element("rma.mv", x$class.f)) 
            res.table <- res.table[, seq_len(8)]
        print(res.table, quote = FALSE, right = TRUE)
    }
    invisible()
}
