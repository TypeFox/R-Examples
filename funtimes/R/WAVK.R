WAVK <- function(z, kn = NULL) {
    if (!is.numeric(z) | !is.vector(z)) {
        stop("input object should be a vector.")
    }    
    if (any(is.na(z))) {
        stop("input vector should not contain missing values.")
    }
    kn <- round(kn)
    T <- length(z)
    ave_group <- sapply(c(1:(T-kn+1)), function(x) mean(z[x:(x+kn-1)]))
    ave_all <- mean(ave_group)
    MST <- sum((ave_group - ave_all)^2) * kn / (T - 1)
    MSE <- sum(sapply(c(1:(T-kn+1)), function(x) sum((z[x:(x+kn-1)] - ave_group[x])^2))) / (T*(kn-1))
    Tn <- MST - MSE
    sigma2 <- sum(diff(z)^2)/(2 * (T - 1))
    Tns <- sqrt(T/kn) * Tn / (sqrt(4/3) * sigma2)
    crit <- pnorm(Tns, mean = 0, sd = 1)
    if (crit < 0.5) {
        p.value <- crit * 2
    } else {
        p.value <- (1 - crit) * 2
    }
    list(Tn = Tn, Tns = Tns, p.value = p.value)
}