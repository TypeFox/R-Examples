# crossings model (Goodman, 1972)
# Ref:
#Goodman, L. (1972).  Some multiplicative models for the analysis of cross-classified data.
#In: Proceedings of the Sixth Berkeley Symposium on Mathematical Statistics and Probability,
#Berkeley, CA: University of California Press, pp. 649-696.

crossings <- function(i, j, n) {
    npar <- n - 1
    result <- list()
    for(c in 1:npar) {
        overi <- c >= i
        overj <- c >= j
        result[[c]] <- (overi & !overj) + (overj & !overi)
    }
    result <- matrix(unlist(result), length(i), npar)
    colnames(result) <- paste('C', 1:npar, sep='')
    result
 }

Crossings <- function(...) {
    dots <- list(...)
    if (length(dots) != 2) stop("Crossings() is defined for only two factors")
    if (length(dots[[1]]) != length(dots[[2]]))
    stop("arguments to Crossings() must all have same length")
    dots <- lapply(dots, as.factor)
    n <- nlevels(dots[[1]])
    if (nlevels(dots[[2]]) != n)
        stop("arguments to Crossings() must all have same number of levels")
    result <- crossings(as.numeric(dots[[1]]), as.numeric(dots[[2]]), n)
    rownames(result) <- do.call("paste", c(dots, sep = ""))
    result
}





