normtol.int <- function (x, alpha = 0.05, P = 0.99, side = 1, method = c("HE", "HE2",
	    "WBE", "ELL", "KM", "EXACT", "OCT"), m = 50, log.norm = FALSE) 
{
    if (log.norm) 
        x <- log(x)
    x.bar <- mean(x)
    s <- sd(x)
    n <- length(x)
    method <- match.arg(method)
    K <- invisible(K.factor(n = n, alpha = alpha, P = P, side = side, 
        method = method, m = m))
    lower <- x.bar - s * K
    upper <- x.bar + s * K
    if (log.norm) {
        lower <- exp(lower)
        upper <- exp(upper)
        x.bar <- exp(x.bar)
    }
    if (side == 1) {
        temp <- data.frame(cbind(alpha, P, x.bar, lower, upper))
        colnames(temp) <- c("alpha", "P", "x.bar", "1-sided.lower", 
            "1-sided.upper")
    }
    else {
        temp <- data.frame(cbind(alpha, P, x.bar, lower, upper))
        colnames(temp) <- c("alpha", "P", "x.bar", "2-sided.lower", 
            "2-sided.upper")
    }
    temp
}








