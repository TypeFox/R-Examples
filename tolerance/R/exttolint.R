exttol.int <- function (x, alpha = 0.05, P = 0.99, side = 1,
    dist = c("Weibull", "Gumbel"), ext = c("min", "max"), NR.delta = 1e-08) 
{
    m.x <- abs(max(x))+1000
    temp.ind <- 0
    if (sum(abs(x)>1000) > 0) {
        temp.ind <- 1
        x <- x/m.x
    }
    if (side != 1 && side != 2) {
        stop(paste("Must specify a one-sided or two-sided procedure!", 
            "\n"))
    }
    if (side == 2) {
        alpha <- alpha/2
    P <- (P + 1)/2
    }
    n <- length(x)
    dist <- match.arg(dist)
    ext <- match.arg(ext)
    if (dist == "Weibull") { 
        x <- log(x)
        ext <- "min"
    }
    delta <- sqrt((mean(x^2) - mean(x)^2) * 6/pi^2)
    x.bar <- mean(x)
    temp <- (dist == "Weibull" | (dist == "Gumbel" & ext == "min"))
    xi <- x.bar + digamma(1) * (1 - 2 * temp)
    theta.old <- c(xi, delta)
    diff <- 1
    if (temp == TRUE) {
        while (sum(diff > NR.delta) > 0) {
            f <- sum(x * exp(x/delta))
            f.1 <- -sum(x^2 * exp(x/delta))/(delta^2)
            g <- sum(exp(x/delta))
            g.1 <- -f/(delta^2)
            d <- delta + x.bar - (f/g)
            d.1 <- 1 - (g * f.1 - f * g.1)/(g^2)
            delta.new <- delta - d/d.1
            xi.new <- -delta.new * log(n/sum(exp(x/delta.new)))
            delta.old <- delta
            xi.old <- xi
            delta <- delta.new
            xi <- xi.new
            if (is.na(xi) | is.na(delta) | delta < 0) {
                xi <- theta.old[1]
                delta <- theta.old[2]
                diff <- NR.delta/5
            }
            else diff <- c(abs(delta.new - delta.old), abs(xi.new - 
                xi.old))
        }
    }
    else {
        lam <- 1/delta
        while (sum(diff > NR.delta) > 0) {
            f <- sum(x * exp(-lam * x))
            f.1 <- -sum(x^2 * exp(-lam * x))
            g <- sum(exp(-lam * x))
            g.1 <- -f
            d <- (1/lam) - x.bar + (f/g)
            d.1 <- (f^2/g^2) + (f.1/g) - (1/lam^2)
            lam.new <- lam - (d/d.1)
            xi.new <- -(1/lam.new) * log((1/n) * sum(exp(-lam.new * 
                x)))
            lam.old <- lam
            xi.old <- xi
            delta.old <- 1/lam
            lam <- lam.new
            xi <- xi.new
            delta.new <- 1/lam
            delta <- delta.new
            if (is.na(xi) | is.na(delta) | delta < 0) {
                xi <- theta.old[1]
                delta <- theta.old[2]
                lam <- 1/delta
                diff <- NR.delta/5
            }
            else diff <- c(abs(delta.new - delta.old), abs(xi.new - 
                xi.old))
        }
    }
    lambda <- function(P) log(-log(P))
    k.t <- function(x1, x2, n) suppressWarnings(qt(1 - x1, df = (n - 
        1), ncp = (-sqrt(n) * lambda(x2))))
    lower <- xi - delta * k.t(alpha, P, n)/sqrt(n - 1)
    upper <- xi - delta * k.t(1 - alpha, 1 - P, n)/sqrt(n - 1)
    if (dist == "Gumbel" & ext == "max") {
#        lower <- xi + delta * k.t(alpha, 1 - P, n)/sqrt(n - 1)
#        upper <- xi + delta * k.t(1 - alpha, P, n)/sqrt(n - 1)
        lower <- xi + delta * k.t(1 - alpha, 1 - P, n)/sqrt(n - 1)
        upper <- xi + delta * k.t(alpha, P, n)/sqrt(n - 1)
    }
    a <- xi
    b <- delta
    if (dist == "Weibull") {
        a <- 1/delta
        b <- exp(xi)
        lower <- exp(lower)
        upper <- exp(upper)
    }
    if (side == 2) {
        alpha <- 2 * alpha
    P <- (2 * P) - 1
    }
    if (temp.ind==1) {
        b <- b*m.x
        lower <- lower*m.x
        upper <- upper*m.x
    }
    temp <- data.frame(cbind(alpha, P, a, b, lower, upper))
    if (side == 2) {
        colnames(temp) <- c("alpha", "P", "shape.1", "shape.2", 
        "2-sided.lower", "2-sided.upper")
    }
    else {
        colnames(temp) <- c("alpha", "P", "shape.1", "shape.2", 
        "1-sided.lower", "1-sided.upper")
    }
    temp
}
