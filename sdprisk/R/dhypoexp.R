dhypoexp <- function(x, rate = 1.0, log = FALSE) {
    stopifnot(is.numeric(x),
              all(is.finite(rate)),
              is.logical(log))

    res <- drop(tcrossprod(outer(X   = x,
                                 Y   = rate,
                                 FUN = dexp,
                                 log = FALSE),
                           t(ratetoalpha(rate))))

    if (log) {
        return(log(res))
    } else {
        return(res)
    }
}
