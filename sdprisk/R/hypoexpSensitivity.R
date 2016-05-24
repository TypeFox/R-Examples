hypoexpSensitivity <- function(process) {
    psi.1 <- hypoexpRuinprob(process)[['psi.1']]
    varu  <- hypoexpVaru(process)
    tvaru <- hypoexpTvaru(process)

    varu.sens <- function(prob) {
        return(-1.0 / (process[['p']] * process[['zeta']] * psi.1(varu(prob))))
    }

    tvaru.sens <- function(prob) {
        return((varu(prob) - tvaru(prob)) / prob)
    }

    return(list(varu.sens  = varu.sens,
                tvaru.sens = tvaru.sens))
}
