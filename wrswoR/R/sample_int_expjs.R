#' @rdname sample_int
#' @examples
#' ## Algorithm A-ExpJ (paper version)
#' s <- sample_int_expjs(20000, 10000, runif(20000))
#' stopifnot(unique(s) == s)
#' p <- c(995, rep(1, 5))
#' n <- 1000
#' set.seed(42)
#' tbl <- table(replicate(sample_int_expjs(6, 3, p),
#'                        n = n)) / n
#' stopifnot(abs(tbl - c(1, rep(0.4, 5))) < 0.04)
#'
"sample_int_expjs"
