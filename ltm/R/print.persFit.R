print.persFit <-
function (x, digits = 4, ...) {
    if (!inherits(x, "persFit"))
        stop("Use only with 'persFit' objects.\n")
    cat("\nPerson-Fit Statistics and P-values\n")
    cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n", sep = "")
    Alt <- switch(x$alternative,
        "less" = "Inconsistent response pattern under the estimated model",
        "greater" = "More consistent response pattern than the model predicts",
        "two.sided" = "Either inconsistent or more consistent response pattern under the estimated model")
    cat("\nAlternative:", Alt)
    if (x$simulate.p.value)
        cat("\nMonte Carlo samples:", x$B, "\n\n")
    else
        cat("\n\n")
    out.dat1 <- as.data.frame(round(x$Tobs, digits))
    out <- apply(x$p.value, 2, function (x) {
        val <- round(x, digits)
        res <- formatC(val, digits = digits)
        res[val == 0] <- paste("<0.", paste(rep(0, digits - 1), collapse = ""), "1", collapse = "", sep = "")
        res
    })
    out.dat2 <- as.data.frame(if (!is.matrix(out)) rbind(out) else out)
    names(out.dat2) <- switch(x$alternative,
        "less" = paste("Pr(<", names(out.dat2), ")", sep = ""),
        "greater" = paste("Pr(>", names(out.dat2), ")", sep = ""),
        "two.sided" = paste("Pr(>|", names(out.dat2), "|)", sep = ""))
    out.dat <- cbind(out.dat1, out.dat2)
    if (length(out.dat) == 4)
        out.dat <- out.dat[c(1, 3, 2, 4)]
    print(cbind(as.data.frame(x$resp.patterns), out.dat))
    cat("\n\n")
    invisible(x)    
}
