exp2tol.int<-function (x, alpha = 0.05, P = 0.99, side = 1, method = c("GPU",
    "DUN","KM"), type.2 = FALSE)
{
    if (side != 1 && side != 2) {
        stop(paste("Must specify a one-sided or two-sided procedure!",
            "\n"))
    }
    method <- match.arg(method)
    if (side == 2) {
        alpha <- alpha/2
        P <- (P + 1)/2
    }
    n <- length(x)
    T <- min(x)
    S <- sum(x - T)
    if (type.2) {
        mx <- max(x)
        r <- n - sum(x == mx)
        m <- function(P, R, n) (1 + n * log(P))/(r - (5/2))
        k1 <- (-m(P, r, n) - qnorm(1 - alpha) * sqrt(m(P, r, 
            n)^2/r + (1/r^2)))/n
        k2 <- (-m(1 - P, r, n) - qnorm(alpha) * sqrt(m(1 - P, 
            r, n)^2/r + (1/r^2)))/n
    }
    else {
        k1 <- (1 - ((P^n)/alpha)^(1/(n - 1)))/n
   if (method == "KM") {
  k2 <- (1 - (((1-P)^n)/(1-alpha))^(1/(n-1)))/n
   } else {
         k2 <- qchisq(P, df = 2)/qchisq(alpha, df = 2 * n - 2)
         if (method == "DUN") {
             lambda = 1.71 + 1.57 * log(log(1/alpha))
             k2 <- k2 - (lambda/n)^(1.63 + 0.39 * lambda)
         }
   }
    }
    lower <- T + S * k1
    upper <- T + S * k2
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