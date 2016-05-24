mean.claiminfo <- function(x, ...) {
    tryCatch(expr  = x[['mu']],
             error = function(.err) NA_real_)
}
