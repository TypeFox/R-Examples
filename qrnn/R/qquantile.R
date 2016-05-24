qquantile <-
function(p, tau, quant, lower = -Inf, tol = .Machine$double.eps^0.25,
         maxiter = 1000, range.mult = 1.1, max.error = 100, ...)
{
    if (length(tau) != length(quant)) 
        stop("\"tau\" and \"quant\" must be same length")
    if (any(tau < 0) | any(tau > 1))
        stop("\"tau\" must be in range [0, 1]")
    if (any(p < 0) | any(p > 1))
        stop("\"p\" must be in range [0, 1]")
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
    q <- approx(x = tau, y = quant, xout = p)$y
    qq.lu <- function(p, tau, quant, min.max, tol, maxiter, ...) {
        cost <- function(q, p, tau, quant, ...){
            pp <- pquantile(q, tau = tau, quant = quant, ...)
            pp-p
        }
        if(p < tau[1]){
            q <- uniroot(f = cost, lower = min(min.max),
                         upper = quant[1], tol = tol,
                         maxiter = maxiter, p = p, tau = tau,
                         quant = quant, ...)$root
        } else {
            q <- uniroot(f = cost, lower = quant[length(quant)],
                         upper = max(min.max), tol = tol,
                         maxiter = maxiter, p = p, tau = tau,
                         quant = quant, ...)$root        
        }
        q
    }
    if (any(is.na(q))) {
        error <- TRUE
        n.error <- 0
        quant.range <- diff(range(quant))
        while (error & (n.error < max.error)) {
            quant.range <- range.mult*quant.range
            min.max <- c(quant[1] - quant.range,
                         quant[length(quant)] + quant.range)
            qq <- try(sapply(p[is.na(q)], qq.lu, tau = tau,
                             quant = quant, min.max = min.max,
                             tol = tol, maxiter = maxiter, ...),
                      silent = TRUE)
            error <- inherits(qq, "try-error")
            n.error <- n.error + 1
        }
        if (error) stop("n.error > max.error")
        q[is.na(q)] <- qq
    }
    if (lower == -Inf) q[p == 0] <- -Inf
    q[p == 1] <- Inf
    q
}
