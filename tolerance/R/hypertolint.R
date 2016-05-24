hypertol.int <- function (x, n, N, m = NULL, alpha = 0.05, P = 0.99, side = 1, 
    method = c("EX", "LS", "CC")) 
{
    if (side != 1 && side != 2) {
        stop(paste("Must specify a one-sided or two-sided procedure!", 
            "\n"))
    }
    if (side == 2) {
        alpha <- alpha/2
        P <- (P + 1)/2
    }
    rate <- n/N
    if (rate < 0.05) 
        warning("Sampling rate < 0.05.  Results may not be accurate!", 
            call. = FALSE)
    method <- match.arg(method)
    if (length(x) > 1) 
        x <- sum(x)
    p.hat <- x/n
    k <- qnorm(1 - alpha)
    if (is.null(m)) 
        m <- n
    fpc <- sqrt((N - n)/(N - 1))
    if (method == "EX") {
				temp <- x:N
				lower.ex <- 1-phyper(x-1,temp,N-temp,n)
				Ml <- temp[min(length(temp),which(lower.ex>alpha))]
				upper.ex <- phyper(x,temp,N-temp,n)
				Mu <- temp[max(1,which(upper.ex>alpha))]
    }
    if (method == "LS" |method == "CC") {
        lower.p <- p.hat - k * sqrt(p.hat * (1 - p.hat)/n) * 
            fpc - 1/(2 * n)*(method == "CC")
        upper.p <- p.hat + k * sqrt(p.hat * (1 - p.hat)/n) * 
            fpc + 1/(2 * n)*(method == "CC")
				lower.p <- max(0, lower.p)
				upper.p <- min(upper.p, 1)
				Ml <- max(0, floor(N * lower.p))
				Mu <- min(ceiling(N * upper.p), N)
    }
    lower <- qhyper(1 - P, m = Ml, n = N - Ml, k = m)
    upper <- qhyper(P, m = Mu, n = N - Mu, k = m)
    if (side == 2) {
        alpha <- 2 * alpha
        P <- (2 * P) - 1
    }
    temp <- data.frame(cbind(alpha, P, rate, p.hat, lower, upper))
    if (side == 2) {
        colnames(temp) <- c("alpha", "P", "rate", "p.hat", "2-sided.lower", 
            "2-sided.upper")
    }
    else {
        colnames(temp) <- c("alpha", "P", "rate", "p.hat", "1-sided.lower", 
            "1-sided.upper")
    }
    temp
}