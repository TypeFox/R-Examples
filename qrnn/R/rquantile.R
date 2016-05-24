rquantile <-
function(n, tau, quant, lower = -Inf, tol = .Machine$double.eps^0.25,
         maxiter = 1000, range.mult = 1.1, max.error = 100, ...)
{
    if (length(tau) != length(quant)) 
        stop("\"tau\" and \"quant\" must be same length")
    if (any(tau < 0) | any(tau > 1))
        stop("\"tau\" must be in range [0, 1]")
    quant[quant < lower] <- lower
    if (is.unsorted(tau) | is.unsorted(quant)) {
        warning("sorting \"tau\" or \"quant\"")
        tau <- sort(tau)
        quant <- sort(quant)
    }
    q <- qquantile(runif(n), tau = tau, quant = quant,
                   lower = lower, tol = tol, maxiter = maxiter,
                   range.mult = range.mult, max.error = max.error,
                   ...)
    q
}
