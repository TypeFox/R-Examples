# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

#' @export
print.cca <- function(x, ...) {
    # print function call
    if(!is.null(call <- x$call)) {
        cat("\nCall:\n")
        dput(x$call)
    }
    # print maximum correlation
    cat("\nCanonical correlations:\n")
    print(x$cor, ...)
    # return object invisibly
    invisible(x)
}

#' @export
print.maxCor <- function(x, ...) {
    # print function call
    if(!is.null(call <- x$call)) {
        cat("\nCall:\n")
        dput(x$call)
    }
    # print maximum association
    cat("\nMaximum correlation:\n")
    print(x$cor, ...)
    # return object invisibly
    invisible(x)
}

#' @export
print.permTest <- function(x, ...) {
    # print general statement
    cat("\nPermutation test for no association\n\n")
    # print maximum correlation and p-value
    cat(sprintf("r = %f, p-value = %f\n", x$cor0, x$pValue))
    # print number of random permuations
    cat(sprintf("R = %d random permuations\n", x$R))
    # print alternative hypothesis
    cat("Alternative hypothesis: true maximum correlation is not equal to 0\n")
    # return object invisibly
    invisible(x)
}
