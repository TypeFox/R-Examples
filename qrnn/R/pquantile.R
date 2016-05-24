pquantile <-
function(q, tau, quant, lower = -Inf, ...)
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
    if ((lower != -Inf) & (!(0 %in% tau))) {
        quant <- c(lower, quant)
        tau <- c(0, tau)
    } else if (lower == -Inf) {
        tau <- tau[quant != -Inf]
        quant <- quant[quant != -Inf]
    }
    p <- q
    p[q >= lower] <- approx(x = quant, y = tau, xout = q[q >= lower],
                           ties = max)$y
    pq.lu <- function(q, tau, quant, ...) {
        if (q < quant[1]){
            p <- integrate(dquantile, lower = -Inf, upper = q,
                           tau = tau, quant = quant, ...)$value
        } else {
            p <- integrate(dquantile, lower = quant[length(quant)],
                           upper = q, tau = tau, quant = quant,
                           ...)$value + tau[length(quant)]            
        }
        p
    }
    if (any(is.na(p))) {
        p[is.na(p)] <- sapply(q[is.na(p)], pq.lu, tau = tau,
                              quant = quant, ...)
    }
    if (any(q < lower)) {
        warning("\"q\" < lower limit")
        p[q < lower] <- NA
    }
    if (lower == -Inf) p[q == -Inf] <- 0
    p[q == Inf] <- 1
    p
}
