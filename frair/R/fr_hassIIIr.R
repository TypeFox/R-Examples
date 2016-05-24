# Hassell's Type III pred-prey function not assuming replacement.
# Now deprecated (gracefully). 
# TODO: Remove in future versions (> 1.0).
hassIIIr <- function(...) {
    warning("hassIIIr has been renamed and deprecated. Use hassIIInr instead.")
    hassIIInr(...)
}

hassIIIr_fit <- function(...) {
    warning("hassIIIr_fit has been renamed and deprecated. Use hassIIInr_fit instead.")
    hassIIInr_fit(...)
}

hassIIIr_nll <- function(...) {
    warning("hassIIIr_nll has been renamed and deprecated. Use hassIIInr_nll instead.")
    hassIIInr_nll(...)
}

hassIIIr_diff <- function(...) {
    warning("hassIIIr_diff has been renamed and deprecated. Use hassIIInr_diff instead.")
    hassIIInr_diff(...)
}

# The NLL for the difference model... used by frair_compare()
hassIIIr_nll_diff <- function(...) {
    warning("hassIIIr_nll_diff has been renamed and deprecated. Use hassIIInr_nll_diff instead.")
    hassIIIr_nll_diff(...)
}
