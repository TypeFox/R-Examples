lundbergRuinprob <- function(process, use.factor = FALSE) {
    stopifnot(is.logical(use.factor))

    myadj <- adjcoef(process)
    if (is.null(myadj)) {
        stop('Unable to compute the adjustment coefficient.')
    }

    const <- 1.0
    if (use.factor) {
        .pdf <- tryCatch(expr  = process$claims$pdf,
                         error = function(.error) NULL)

        if (!is.null(.pdf)) {
            const <- process[['p']] / (myadj * (process[['q']] / mean(process[['claims']]) *
                      integrate(f = function(.x) {
                                    .res <- .pdf(.x) * .x * exp(.x * myadj)
                                    .res[which(is.nan(.res))] <- 0.0
                                    return(.res)
                                },
                                lower = 0.0,
                                upper = Inf)$value
                      + 1.0 / process[['zeta']]
                      ))
        } else {
            warning('Cannot compute the coefficient without the claim density.\n',
                    'Using 1 as a replacement as with use.factor = FALSE.')
        }
    }

    return(structure(.Data = function(x) const * exp(-x * myadj),
                     C     = const,
                     r     = myadj))
}
