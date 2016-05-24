# 5) c constant for Tukey Biweight S
Tbsc <- function(alpha, q, eps=10^-8, maxit=1000) {
    ctest <- talpha <- sqrt(qchisq(1 - alpha, q))
    diff <- Inf
    iter <- 1
    while ((diff > eps) && (iter < maxit)) {
        cold <- ctest
        ctest <- Tbsb(cold, q)/alpha
        diff <- abs(cold - ctest)
        iter <- iter + 1
    }
    ctest
}
