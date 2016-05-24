varu <- function(process, method = c('saddlepoint', 'hypoexp'), ...) {
    stopifnot(is.riskproc(process))
    if (missing(method)) {
        if (is.hypoexp(process)) {
            method <- 'hypoexp'
        } else {
            method <- 'saddlepoint'
        }
    } else {
        method <- match.arg(method)
    }

    switch(method,
        saddlepoint = saddlepointVaru(process, ...),
        hypoexp     = hypoexpVaru(process, ...)
    )
}
