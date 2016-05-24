#' @S3method print treemox.pls
print.treemox.pls <-
function(x, ...)
{
    cat("PLS-PM results of Segmentation Trees ", "\n\n")
    cat("$weights", "\n")
    print(x$weights, print.gap=2)
    cat("\n")
    cat("$loadings", "\n")
    print(x$loadings, print.gap=2)
    cat("\n")
    cat("$paths", "\n")
    print(x$paths, print.gap=2)
    cat("\n")
    cat("$r2", "\n")
    print(x$r2, print.gap=2)
    cat("\n")
    invisible(x)
}

