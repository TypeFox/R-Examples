effectBias <- function(y, u, w, theta, eta){

    n <- length(y)
    k <- 1 + floor(n / 2)
    ratio <- abs(y) / u
    p <- 2 * pnorm(-ratio)

    bias <- rep(NA, n)
    toIntBias <- function(z){return(z * dnorm(z))}
    b0 <- 0

    for (i in 1:n){
        Bj <- rep(NA, k)
        Cj <- rep(NA, k)
        b <- -u[i] * qnorm(p / 2)

        # j = 1
        j <- 1
        I1 <- integrate(toIntBias, lower = (-b[2 * j] - theta) / eta[i], upper = (b[2 * j] - theta) / eta[i])$value
        I2 <- integrate(toIntBias, lower = (b0 - theta) / eta[i], upper = (-b0 - theta) / eta[i])$value
        Bj[j] <- I1 + I2
        Cj[j] <- pnorm((b[2 * j] - theta) / eta[i]) - pnorm((-b[2 * j] - theta) / eta[i]) + pnorm((-b0 - theta) / eta[i]) - pnorm((b0 - theta) / eta[i])

        # j = 2, ..., k - 1
        for (j in 2:(k - 1)){
            I1 <- integrate(toIntBias, lower = (-b[2 * j] - theta) / eta[i], upper = (b[2 * j] - theta) / eta[i])$value
            I2 <- integrate(toIntBias, lower = (b[2 * j - 2] - theta) / eta[i], upper = (-b[2 * j - 2] - theta) / eta[i])$value
            Bj[j] <- I1 + I2
            Cj[j] <- pnorm((b[2 * j] - theta) / eta[i]) - pnorm((-b[2 * j] - theta) / eta[i]) + pnorm((-b[2 * j - 2] - theta) / eta[i]) - pnorm((b[2 * j - 2] - theta) / eta[i])
        }

        # j = k
        j <- k
        I1 <- integrate(toIntBias, lower = -Inf, upper = Inf)$value
        I2 <- integrate(toIntBias, lower = (b[2 * j - 2] - theta) / eta[i], upper = (-b[2 * j - 2] - theta) / eta[i])$value
        Bj[j] <- I1 + I2
        Cj[j] <- 1 + pnorm((-b[2 * j - 2] - theta) / eta[i]) - pnorm((b[2 * j - 2] - theta) / eta[i])

        # compute expected bias using estimated w
        bias[i] <- eta[i] * sum(w * Bj) / sum(w * Cj)
    }
    
    ## compute bias and generate output
    bias <- sign(y) * bias
    dat <- round(cbind(y, u, p, bias, y - bias, bias / y), 4)
    colnames(dat) <- c("y", "u", "p", "bias", "y - bias", "bias / y")
    
    return(list("dat" = dat))
}
