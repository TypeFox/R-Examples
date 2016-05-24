uniftol.int <- function (x, alpha = 0.05, P = 0.99, upper = NULL, lower = NULL,
    side = 1) 
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
    x.1 <- lower
    x.n <- upper
    if (is.null(x.1)) 
        x.1 <- min(x)
    if (is.null(x.n)) 
        x.n <- max(x)
    lower <- ((1 - P)/(1 - alpha)^(1/n)) * (x.n - x.1) + x.1
    upper <- (P/(alpha)^(1/n)) * (x.n - x.1) + x.1
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
