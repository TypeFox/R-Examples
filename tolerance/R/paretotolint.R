paretotol.int <- function (x, alpha = 0.05, P = 0.99, side = 1, method = c("GPU", 
    "DUN"), power.dist = FALSE) 
{
    if (side != 1 && side != 2) {
        stop(paste("Must specify a one-sided or two-sided procedure!", 
            "\n"))
    }
    if (power.dist) {
        x <- log(1/x)
    }
    else {
        x <- log(x)
    }
    method <- match.arg(method)
    out <- exp2tol.int(x = x, alpha = alpha, P = P, side = side, method = method,
            type.2 = FALSE)
    lower <- out[1, 3]
    upper <- out[1, 4]
    if (power.dist) {
        lower.1 <- 1/exp(upper)
        upper <- 1/exp(lower)
        lower <- lower.1
    }
    else {
        lower <- exp(lower)
        upper <- exp(upper)
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
