dquantile <-
function (x, tau, quant, lower = -Inf) 
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
    if (any(x < lower)) {
        warning("\"x\" < lower limit; replacing values with \"lower\"")
        x[x < lower] <- lower
    }
    if ((lower != -Inf) & (!(0 %in% tau))) {
        quant <- c(lower, quant)
        tau <- c(0, tau)
    }
    dq <- function(x, tau, quant) {
        n <- length(quant)
        d <- 0
        z1 <- (tau[2]-tau[1])/(quant[2]-quant[1])
        b1 <- tau[1]/z1
        zn <- (tau[n]-tau[n-1])/(quant[n]-quant[n-1])
        bn <- (1-tau[n])/zn
        run <- TRUE
        if (x < quant[1]) {
            d <- z1*exp(-abs(x-quant[1])/b1)
            run <- FALSE
        }
        if (x >= quant[n]) {
            d <- zn*exp(-abs(x-quant[n])/bn)
            run <- FALSE
        }
        if (run) {
            for (j in 2:n) {
                if (quant[j] > x) {
                  d <- (tau[j]-tau[j-1])/(quant[j]-quant[j-1])
                  break
                }
            }
        }
        d
    }
    d <- sapply(x, dq, tau = tau, quant = quant)
    if((lower != -Inf) & any(x == lower)) {
        d[x == lower] <- max(tau[quant == lower])
    }
    d
}
