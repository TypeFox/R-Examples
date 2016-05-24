print.pred.density <-
function (x, digits = NULL, ...) 
{
    outmat = matrix(numeric(0), length(x$fit), 2)
    colnames(outmat) = c("Exp.Val.", "Std.Err.")
    rownames(outmat) = names(x$fit)
    outmat[, 1] = x$fit
    outmat[, 2] = x$std.err
    cat("Call:\n")
    print(x$call)
    cat(paste("\nDensities for conditional forecast(s)\n", x$n, 
        " data points, based on ", x$nmodel, " models;\n", sep = ""))
    print(outmat, digits = digits, ...)
}
