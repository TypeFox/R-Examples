pvaldistance <- function(x, method = c('ks', 'cvm'), dist.to = c('uniform')) {
    stopifnot(is.numeric(x), is.character(method), is.character(dist.to))

    method  <- match.arg(method)
    dist.to <- match.arg(dist.to)

    if (dist.to == 'uniform') {
        x   <- sort(as.vector(x))
        num <- length(x)

        switch(method,
            ks  = max(abs(sweep(x      = cbind(c(0.0, x), c(x, 1.0)),
                                MARGIN = 1L,
                                STATS  = seq.int(0.0, 1.0, 1.0 / num)))),
            cvm = sum((seq_len(num) / num - 0.5 / num - x)^2.0) + 1.0 / (12.0 * num)
        )
    } else {
        NA_real_
    }
}
