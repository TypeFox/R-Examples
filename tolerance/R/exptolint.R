exptol.int <- function (x, alpha = 0.05, P = 0.99, side = 1, type.2 = FALSE) 
{
    if (side != 1 && side != 2) {
        stop(paste("Must specify a one-sided or two-sided procedure!", 
            "\n"))
    }
    if (side == 2) 
        alpha <- alpha/2
    n <- length(x)
    l.hat <- mean(x)
    if (type.2) {
        mx <- max(x)
        r <- n - sum(x == mx)
    }
    else r <- n
    if (side == 2) {
        lower <- 2 * r * l.hat * log(2/(1 + P))/qchisq(1 - alpha, 
            df = 2 * r)
        upper <- 2 * r * l.hat * log(2/(1 - P))/qchisq(alpha, 
            df = 2 * r)
        alpha <- 2 * alpha
    }
    else {
        lower <- 2 * r * l.hat * log(1/P)/qchisq(1 - alpha, df = 2 * 
            r)
        upper <- 2 * r * l.hat * log(1/(1 - P))/qchisq(alpha, 
            df = 2 * r)
    }
    temp <- data.frame(cbind(alpha, P, l.hat, lower, upper))
    if (side == 2) {
        colnames(temp) <- c("alpha", "P", "lambda.hat", "2-sided.lower", 
            "2-sided.upper")
    }
    else {
        colnames(temp) <- c("alpha", "P", "lambda.hat", "1-sided.lower", 
            "1-sided.upper")
    }
    temp
}
