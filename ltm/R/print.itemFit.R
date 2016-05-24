print.itemFit <-
function (x, digits = 4, ...) {
    if (!inherits(x, "itemFit"))
        stop("Use only with 'itemFit' objects.\n")
    cat("\nItem-Fit Statistics and P-values\n")
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
    cat("\nAlternative: Items do not fit the model")
    cat("\nAbility Categories:", x$G)
    if (x$simulate.p.value)
        cat("\nMonte Carlo samples:", x$B, "\n\n")
    else
        cat("\n\n")
    pvals <- round(x$p.values, digits = digits)
    out.pvals <- formatC(pvals, digits = digits)
    out.pvals[pvals == 0] <- paste("<0.", paste(rep(0, digits - 1), collapse = ""), "1", collapse = "", sep = "")
    print(data.frame("X^2" = round(x$Tobs, digits), "Pr(>X^2)" = out.pvals, row.names = names(x$Tobs), 
            check.names = FALSE))
    cat("\n\n")
    invisible(x)    
}
