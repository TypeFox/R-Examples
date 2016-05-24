#' @export
pmfhistogram <-
function (dist, x = 0:10, argsForDist = list(), breaks, scale, 
    ...) 
{
    gotBreaks <- !missing(breaks)
    if (is.integer(x) && !gotBreaks) {
        breaks = seq(min(x) - 0.5, max(x) + 0.5, by = 1)
        gotBreaks <- TRUE
    }
    probs <- do.call(dist, args = c(list(x = x), argsForDist))
    if (missing(scale)) {
        scale = min(50000, 5/min(probs))
    }
    x = rep(x, times = round(probs * scale))
    if (gotBreaks) {
        histogram(~x, breaks = breaks, ...)
    }
    else {
        histogram(~x, ...)
    }
}
