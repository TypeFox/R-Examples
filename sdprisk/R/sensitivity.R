sensitivity <- function(process, method = c('saddlepoint', 'hypoexp'), ...) {
    stopifnot(is.riskproc(process))
    switch(match.arg(method),
           saddlepoint = saddlepointSensitivity(process, ...),
           hypoexp     = hypoexpSensitivity(process, ...))
}
