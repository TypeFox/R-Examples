#solve.Rcpp_Rmumps <- function(a,b) if(missing(b)) a$inv() else a$solve(b)
solve.Rcpp_Rmumps <- function (a, b, ...) if(missing(b)) a$inv() else a$solve(b)
dim.Rcpp_Rmumps <- function (x) x$dim()
nrow.Rcpp_Rmumps <- function (x) x$nrow()
ncol.Rcpp_Rmumps <- function (x) x$ncol()
print.Rcpp_Rmumps <- function (x, ...) x$print()
show.Rcpp_Rmumps <- function (x) x$show()
#setMethod("solve",
#   signature(a = "Rcpp_Rmumps", b="ANY"),
#   solve.Rcpp_Rmumps
#)
#setMethod("solve",
#   signature(a = "Rcpp_Rmumps", b="Matrix"),
#   solve.Rcpp_Rmumps
#)
#setMethod("solve",
#   signature(a = "Rcpp_Rmumps", b="simple_triplet_matrix"),
#   solve.Rcpp_Rmumps
#)
