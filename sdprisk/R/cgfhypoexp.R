cgfhypoexp <- function(x, rate = 1.0, difforder = 0L) {
    difforder <- match(difforder, 0L:2L, 4L)
    is.na(x)  <- (x >= min(rate))
    m <- Map(f         = mgfhypoexp,
             difforder = 0L:2L,
             MoreArgs  = list(x = x, rate = rate))
    switch(difforder,
           log(m[[1L]]),
           m[[2L]] / m[[1L]],
           (m[[3L]] * m[[1L]] - m[[2L]]^2.0) / m[[1L]]^2.0,
           rep.int(NA_real_, length(x)))
}
