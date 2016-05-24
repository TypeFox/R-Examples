hypoexpVaru <- function(process) {
    ruinprob <- hypoexpRuinprob(process)
    psi.r    <- attr(ruinprob, 'diagnostics')[['r']]
    psi      <- ruinprob[['psi']]

    if (length(psi.r) == 1L) {
        varu <- function(prob) {
            (log(attr(ruinprob, 'diagnostics')[['C']]) - log(prob)) / psi.r
        }
    } else {
        varu <- function(prob) {
            ## Compute the Lundberg approximation (without multiplicative constant),
            ## multiplied by 4.0 (more generally: a constant bigger than 1) to err on the safe
            ## side, and use this a upper bound for the interval to be considered for the
            ## VaRu.

            vapply(X         = mapply(FUN      = uniroot,
                                      upper    = -4.0 * log(prob) / adjcoef(process),
                                      .prob    = prob,
                                      MoreArgs = list(f     = function(x, .prob) psi(x) - .prob,
                                                      lower = 0.0),
                                      SIMPLIFY = FALSE),
                   FUN       = `[[`,
                   FUN.VALUE = numeric(1L),
                   index     = 'root')
        }
    }

    return(varu)
}
