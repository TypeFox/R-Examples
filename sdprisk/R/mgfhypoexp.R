mgfhypoexp <- function(x, rate = 1.0, difforder = 0L) {
    difforder <- match(as.integer(difforder), 0L:2L, 4L)
    is.na(x)  <- (x >= min(rate))
    inv.diff  <- 1.0 / outer(rate, x, '-')
    prod(rate) * apply(inv.diff, 2L, prod) * switch(difforder,
        1.0,
        colSums(inv.diff),
        colSums(inv.diff)^2.0 + colSums(inv.diff^2.0),
        rep.int(NA_real_, length(x))
    )
}
