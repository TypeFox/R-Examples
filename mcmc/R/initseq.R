initseq <- function(x) {
    stopifnot(is.numeric(x))
    stopifnot(is.finite(x))
    .Call("initseq", x - mean(x), PACKAGE = "mcmc")
}
