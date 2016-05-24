## Wrapper function for easier simulation
rpsim <- function(model) {
    stopifnot(names(model) %in% c('premium', 'freq', 'variance', 'rates'))
    function(timehoriz, num = 8L, m = 1001L) {
        num <- as.integer(num)
        rriskproc(m        = m,
                  window   = c(0.0, timehoriz),
                  num      = num,
                  sigma    = sqrt(model$variance),
                  freq     = model$freq,
                  drift    = model$premium,
                  jumpdist = rhypoexp,
                  rate     = model$rates)
    }
}
