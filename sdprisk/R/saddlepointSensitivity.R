saddlepointSensitivity <- function(process, ...) {

    varu  <- saddlepointVaru(process)
    tvaru <- saddlepointTvaru(process, ...)

    const <- attr(saddlepointRuinprob(process = process, normalize = TRUE, jensen = FALSE), 'diagnostics')[['corrconst']]

    return(list(varu.sens  = function(prob, normalize = TRUE) {
                    v <- attr(varu(prob), 'saddlepoint')
                    return(-ifelse(normalize, const, 1.0) * sqrt(process[['KL.d2']](v)) / dnorm(process[['rv']](v)))
                },
                tvaru.sens = function(prob) {
                    return((varu(prob) - tvaru(prob)) / prob)
                }))
}
