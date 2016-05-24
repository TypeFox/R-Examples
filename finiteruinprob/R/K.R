## CGF (and its derivatives) of the aggregate loss process at time 1.0
K <- function(model) {
    function(x, difforder = 0L) {
        mx <- sdprisk::mgfhypoexp(x         = x,
                                  rate      = model$rates,
                                  difforder = difforder)

        difforder <- match(difforder, 0L:2L, 4L)
        is.na(x)  <- vapply(x >= min(model$rates), isTRUE, logical(1L))

        drop(crossprod(c(model$freq, model$premium, model$variance),
                       switch(difforder,
                              rbind(mx - 1.0,    -x,   0.5 * x^2.0),
                              rbind(mx,        -1.0,             x),
                              rbind(mx,         0.0,           1.0))))
    }
}
