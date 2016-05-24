laptol.int <- function (x, alpha = 0.05, P = 0.99, side = 1) 
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
    k.b <- log(2 * (1 - P))
    mu.hat <- median(x)
    beta.hat <- mean(abs(x - median(x)))
    k <- (-n * k.b + qnorm(1 - alpha) * sqrt(n * (1 + k.b^2) - 
        qnorm(1 - alpha)^2))/(n - qnorm(1 - alpha)^2)
    lower <- mu.hat - k * beta.hat
    upper <- mu.hat + k * beta.hat
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
