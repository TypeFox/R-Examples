ruinprob.finite.dsim <- function(Z) {
    stopifnot(is(Z, 'mts'))

    Y <- apply(Z, 2L, cummax) - Z

    function(x) {
        structure(sweep(apply(rowMeans(outer(X   = x,
                                             Y   = Y,
                                             FUN = '<='),
                                       na.rm = TRUE,
                                       dims  = 2L),
                              MARGIN = 1L,
                              FUN    = cumsum),
                        MARGIN = 1L,
                        STATS  = time(Z),
                        FUN    = '/') * deltat(Z),
                  dimnames = list(time(Y), x))
    }
}
