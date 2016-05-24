ratetoalpha <- function(rate) {
    stopifnot(is.numeric(rate), all(is.finite(rate)))

    vapply(X   = seq_along(rate),
           FUN = function(i) {
               prod(rate[-i] / (rate[-i] - rate[i]))
           },
           FUN.VALUE = numeric(1L))
}
