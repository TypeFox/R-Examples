adjcoef <- function(process) {
    stopifnot(is.riskproc(process))
    tryCatch(expr  = process[['adjcoef']],
             error = function(e) NA_real_)
}
