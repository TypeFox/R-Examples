print.permutest.rma.uni <-
function (x, digits, signif.stars = getOption("show.signif.stars"), 
    signif.legend = signif.stars, ...) 
{
    if (!is.element("permutest.rma.uni", class(x))) 
        stop("Argument 'x' must be an object of class \"permutest.rma.uni\".")
    if (missing(digits)) 
        digits <- x$digits
    cat("\n")
    if (!x$int.only) {
        cat("Test of Moderators (coefficient(s) ", paste(x$btt, 
            collapse = ","), "): \n", sep = "")
        if (x$knha) {
            cat("F(df1 = ", x$m, ", df2 = ", x$dfs, ") = ", formatC(x$QM, 
                digits = digits, format = "f"), ", p-val* ", 
                .pval(x$QMp, digits = digits, showeq = TRUE, 
                  sep = " "), "\n\n", sep = "")
        }
        else {
            cat("QM(df = ", x$m, ") = ", formatC(x$QM, digits = digits, 
                format = "f"), ", p-val* ", .pval(x$QMp, digits = digits, 
                showeq = TRUE, sep = " "), "\n\n", sep = "")
        }
    }
    res.table <- cbind(estimate = x$b, se = x$se, zval = x$zval, 
        `pval*` = x$pval, ci.lb = x$ci.lb, ci.ub = x$ci.ub)
    colnames(res.table) <- c("estimate", "se", "zval", "pval*", 
        "ci.lb", "ci.ub")
    if (x$knha) 
        colnames(res.table)[3] <- "tval"
    signif <- symnum(x$pval, corr = FALSE, na = FALSE, cutpoints = c(0, 
        0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", 
        "*", ".", " "))
    if (signif.stars) {
        res.table <- cbind(formatC(res.table, digits = digits, 
            format = "f"), signif)
        colnames(res.table)[7] <- ""
    }
    else {
        res.table <- formatC(res.table, digits = digits, format = "f")
    }
    res.table[, 4] <- .pval(x$pval, digits = digits)
    cat("Model Results:")
    cat("\n\n")
    print(res.table, quote = FALSE, right = TRUE, print.gap = 2)
    cat("\n")
    if (signif.legend) 
        cat("---\nSignif. codes: ", attr(signif, "legend"), "\n\n")
    invisible()
}
