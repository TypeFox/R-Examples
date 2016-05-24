##
##  w h i t t a k e r . R  Whittaker Smoothing
##


whittaker <- function(y, lambda = 1600, d = 2) {
    warning("Not yet fully implemented (because of some strange behavior).
  See the example in 'savgol' for working code requiring the SparseM package.")
}

# whittaker <- function(y, lambda = 1600, d = 2) {
#     stopifnot(is.numeric(y))
# 
#     success <- require("SparseM", warn.conflicts = FALSE, quietly = TRUE)
#     # success <- library("SparseM", pos = "package:base",
#     #                    logical.return = TRUE, warn.conflicts = FALSE)
#     if (!success)
#       stop("Function 'whittaker' requires package 'SparseM' to be installed.")
# 
#     m <- length(y)
#     E <- as(m, "matrix.diag.csr")
#     class(E) <- "matrix.csr"
# 
#     Dmat <- diff(E, differences = d)
#     B <- E + (lambda * t(Dmat) %*% Dmat)
#     z <- solve(B, y)
# 
#     return(z)
# }
