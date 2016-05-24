#' @export
multirandom <-
function (n = 100, samples = 6, rdist = rnorm, ...) 
{
    DOTS <- list(...)
    rdist <- match.fun(rdist)
    RFUN <- function() {
        args = DOTS
        args$n = n
        do.call("rdist", args)
    }
    data.frame(value = c(replicate(samples, RFUN())), sample = rep(1:samples, 
        each = n))
}
