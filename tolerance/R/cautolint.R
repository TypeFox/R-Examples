cautol.int <- function (x, alpha = 0.05, P = 0.99, side = 1) 
{
    if (side != 1 && side != 2) {
        stop(paste("Must specify a one-sided or two-sided procedure!", 
            "\n"))
    }
    if (side == 2) {
        alpha <- alpha/2
	P <- (P + 1)/2
    }
    n <- length(x)
    inits <- c(median(x), IQR(x)/2)
    cau.ll <- function(x, pars) sum(-dcauchy(x, location = pars[1], 
        scale = pars[2], log = TRUE))
    out <- suppressWarnings(nlm(cau.ll, p = inits, x = x)$estimate)
    theta.hat <- out[1]
    sigma.hat <- out[2]
    c.factor <- 2 + 2 * (qcauchy(1 - P))^2
    k <- sqrt(c.factor/n) * qnorm(1 - alpha) - qcauchy(1 - P)
    lower <- theta.hat - k * sigma.hat
    upper <- theta.hat + k * sigma.hat
    if (side == 2) {
        alpha <- 2 * alpha
	P <- (2 * P) - 1
    }
    temp <- data.frame(cbind(alpha, P, lower, upper))
    if (side == 2) {
        colnames(temp) <- c("alpha", "P", "2-sided.lower", "2-sided.upper")
    }
    else {
        colnames(temp) <- c("alpha", "P", "1-sided.lower", "1-sided.upper")
    }
    temp
}
