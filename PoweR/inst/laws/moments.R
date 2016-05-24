moments1 <- function(mu, b) {
## Laplace
    expectation <- mu
    variance <- 2 * b ^ 2
    return(list(expectation = expectation, variance = variance))
}

moments2 <- function(mu, sigma) {
## Normal
    expectation <- mu
    variance <- sigma ^ 2
    return(list(expectation = expectation, variance = variance))
}

moments3 <- function(l, s) {
## Cauchy
    expectation <- NA
    variance <- NA
    return(list(expectation = expectation, variance = variance))
}

moments4 <- function(mu, s) {
## Logistic
    expectation <- mu
    variance <- pi ^ 2 * s ^ 2 / 3
    return(list(expectation = expectation, variance = variance))
}

moments5 <- function(a, b) {
## Gamma
    expectation <- a / b
    variance <- a / b ^ 2
    return(list(expectation = expectation, variance = variance))
}

moments6 <- function(alpha, beta) {
## Beta
    expectation <- alpha / (alpha + beta)
    variance <- alpha * beta / (alpha + beta) ^ 2 / (alpha + beta + 1)
    return(list(expectation = expectation, variance = variance))
}

moments7 <- function(a, b) {
## Uniform
    expectation <- (a + b) / 2
    variance <- (b - a) ^ 2 / 12
    return(list(expectation = expectation, variance = variance))
}

moments8 <- function(k) {
## Student
    expectation <- if (k > 1) 0 else NA 
    variance <- if (k <= 2) Inf else k / (k - 2)
    return(list(expectation = expectation, variance = variance))
}

moments9 <- function(k) {
## Chi-squared
    expectation <- k
    variance <- 2 * k
   return(list(expectation = expectation, variance = variance))
}

moments10 <- function(mu, sigma) {
## Log Normal
    expectation <- exp(mu + sigma ^ 2 / 2)
    variance <- (exp(sigma ^ 2) - 1) * exp(2 * mu + sigma ^ 2)
    return(list(expectation = expectation, variance = variance))
}

moments11 <- function(lambda, k) {
## Weibull
    expectation <- k * gamma(1 + 1 / lambda)
    variance <- k ^ 2 * gamma(1 + 2 / lambda) - expectation ^ 2
    return(list(expectation = expectation, variance = variance))
}

moments12 <- function(l, b) {
## Shifted Exponential
    expectation <- l + 1 / b
    variance <- 1 / b ^ 2
    return(list(expectation = expectation, variance = variance))
}

moments13 <- function(j) {
## Power Uniform
    expectation <- 1 / (j + 2)
    variance <- (j + 1) ^ 2 / (2 * j + 3) / (j + 2) ^ 2
    return(list(expectation = expectation, variance = variance))
}

moments14 <- function(k, a, b) {
## Average Uniform
    expectation <- (a + b) / 2
    variance <- (b - a) ^ 2 / (12 * k)
    return(list(expectation = expectation, variance = variance))
}

moments15 <- function(j) {
## UUniform
    expectation <- 0.5
    variance <- (2 * j ^ 2 + 3 * j + 2) / (2 * (2 * j + 3) * (2 * j + 4))
    return(list(expectation = expectation, variance = variance))
}

moments16 <- function(j) {
## VUniform
    expectation <- 0.5
    somme <- 0
    for (k in 0:(j+1)) {
        signe <- sign(k / (j + 1) - 0.5)
        if (signe == 0) signe <- -1
        somme <- somme + (-1) ^ k * choose(j + 1, k) * {(-1) ^ (j + 1) * k ^ (j + 2) / ((j + 1) * (j + 2)) - signe * ((j + 1) / 2 - k) ^ (j + 1) * (((j + 1) / 2 - k) / (j + 2) + k / (j + 1)) }
  }
    res <- somme
    variance <- 1 / (12 * (j + 1)) - 1 / 4 + res / factorial(j + 1) 
    return(list(expectation = expectation, variance = variance))
}

moments17 <- function(mu, sigma, nu, tau) {
## Johnson SU
    expectation <- mu
    variance <- sigma ^ 2
    return(list(expectation = expectation, variance = variance))
}

moments18 <- function(l) {
## Symmetrical Tukey
    expectation <- 0
    variance <- 2 * (1 / (2 * l + 1) - (gamma(l + 1)) ^ 2 / gamma(2 * l + 2)) / l ^ 2
    return(list(expectation = expectation, variance = variance))
}

moments19 <- function(p, m) {
## Location contaminated
    expectation <- p * m
    variance <- 1 - (p * m) ^ 2 + p * m ^ 2
    return(list(expectation = expectation, variance = variance))
}

moments20 <- function(g, d) {
## Johnson SB
    expectation <- NA
    variance <- NA
    return(list(expectation = expectation, variance = variance))
}

moments21 <- function(xi, omega, alpha) {
## Skew Normal
    delta <- alpha / sqrt(1 + alpha ^ 2)
    expectation <- xi + omega * sqrt(2 / pi) * delta
    variance <- omega ^ 2 * (1 - 2 * delta ^ 2 / pi)
    return(list(expectation = expectation, variance = variance))
}

moments22 <- function(p, d) {
## Scale contaminated
    expectation <- 0
    variance <- p * d ^ 2 + 1 - p
    return(list(expectation = expectation, variance = variance))
}

moments23 <- function(mu, sigma, xi) {
## Generalized Pareto
    expectation <- if (xi < 1) mu + sigma / (1 + xi) else NA
    variance <- if (xi < 0.5) sigma ^ 2 / (1 + xi) ^ 2 / (1 + 2 * xi) else NA
    return(list(expectation = expectation, variance = variance))
}

moments24 <- function(mu, sigma, p) {
## Generalized Error Distribution
    expectation <- mu
    variance <- sigma ^ 2 * gamma(3 / p) / gamma(1 / p)
    return(list(expectation = expectation, variance = variance))
}

moments25 <- function(alpha, beta, c, mu) {
## Stable
    expectation <- if (alpha > 1) mu else NA
    variance <- if (alpha == 2) 2 * c ^ 2 else Inf
    return(list(expectation = expectation, variance = variance))
}

moments26 <- function(mu, sigma) {
## Gumbel
    expectation <- mu - sigma * digamma(1)
    variance <- pi ^ 2 * sigma ^ 2 / 6
    return(list(expectation = expectation, variance = variance))
}

moments27 <- function(mu, sigma, alpha) {
## Fréchet
    expectation <- if (alpha > 1) mu + sigma * gamma(1 - 1 / alpha) else Inf
    variance <- if (alpha > 2) sigma ^ 2 * (gamma(1 - 2 / alpha) - (gamma(1 - 1 / alpha)) ^ 2) else Inf
    return(list(expectation = expectation, variance = variance))
}

moments28 <- function(mu, sigma, xi) {
## Generalized Extreme Value
    expectation <- if (xi != 0 && xi < 1 ) mu + sigma * (gamma(1 - xi) - 1) / xi else if (xi == 0) mu - sigma * digamma(1) else Inf
    g1 <- gamma(1 - xi) ; g2 <- gamma(1 - 2 * xi)
    variance <-  if (xi != 0 && xi < 0.5 ) sigma ^ 2 * (g2 - g1 ^ 2) / xi ^ 2 else if (xi == 0) sigma ^ 2 * pi ^ 2 / 6 else Inf
    return(list(expectation = expectation, variance = variance))
}

moments29 <- function(alpha) {
## Geneeralized Arcsine
    expectation <- 1 - alpha
    variance <- (1 - alpha) * alpha / 2
    return(list(expectation = expectation, variance = variance))
}

moments30 <- function(mu, sigma) {
## Folded Normal
    expectation <- sigma * sqrt(2 / pi) * exp(-mu ^ 2 / (2 * sigma ^ 2)) + mu * (1 - 2 * pnorm(-mu / sigma))
    variance <- mu ^ 2 + sigma ^ 2 - (sigma * sqrt(2 / pi) * exp(-mu ^ 2 / (2 * sigma ^ 2)) + mu * (1 - 2 * pnorm(-mu / sigma))) ^ 2
    return(list(expectation = expectation, variance = variance))
}

moments31 <- function(p, m, d) {
## Mixture Normal
    expectation <- m * p
    variance <- (1 - p) * (1 + p * m ^ 2) + p * d ^ 2
    return(list(expectation = expectation, variance = variance))
}

moments32 <- function(a, b) {
## Truncated Normal
    expectation <- (dnorm(a) - dnorm(b)) / (pnorm(b) - pnorm(a))
    variance <- 1 + (a * dnorm(a) - b * dnorm(b)) / (pnorm(b) - pnorm(a)) - ((dnorm(a) - dnorm(b)) / (pnorm(b) - pnorm(a))) ^ 2
    return(list(expectation = expectation, variance = variance))
}

moments33 <- function(a) {
## Normal with outliers
    expectation <- 0
    variance <- 1
    return(list(expectation = expectation, variance = variance))
}

moments34 <- function(t1, t2, t3) {
## Generalized Exponential Power
    expectation <- NA
    variance <- NA
    return(list(expectation = expectation, variance = variance))
}

moments35 <- function(lambda) {
## Exponential
    expectation <- 1 / lambda
    variance <- 1 / lambda ^ 2
    return(list(expectation = expectation, variance = variance))
}

moments36 <- function(mu, b, k) {
## Asymmetric Laplace
    expectation <- mu + b * (1 / k - k) / sqrt(2)
    variance <- b ^ 2 * (1 + k ^ 4) / (2 * k ^2)
    return(list(expectation = expectation, variance = variance))
}

moments37 <- function(alpha, beta, delta, mu) {
## Normal Inverse Gaussian
    gamma <- sqrt(alpha ^ 2 - beta ^ 2)
    expectation <- mu + beta * delta / gamma
    variance <- delta * alpha ^ 2 / gamma ^ 3
    return(list(expectation = expectation, variance = variance))
}

moments38 <- function(theta, phi, alpha, lambda) {
## Asymmetric Power Distribution
    delta <- 2 * alpha ^ lambda * (1 - alpha) ^ lambda / (alpha ^ lambda + (1 - alpha) ^ lambda)
    expectation <- theta + phi * gamma(2 / lambda) * (1 - 2 * alpha) * delta ^ (-1 / lambda) / gamma(1 / lambda)
    variance <- phi ^ 2 * (gamma(3 / lambda) * gamma(1 / lambda) * (1 - 3 * alpha + 3 * alpha ^ 2) - (gamma(2 / lambda)) ^ 2 * (1 - 2 * alpha) ^ 2) * delta ^ (-2 / lambda) / (gamma(1 / lambda)) ^ 2
    return(list(expectation = expectation, variance = variance))
}

moments39 <- function(mu, sigma, theta1, theta2) {
## modified Asymmetric Power Distribution
    delta <- 2 * theta1 ^ theta2 * (1 - theta1) ^ theta2 / (theta1 ^ theta2 + (1 - theta1) ^ theta2)
    expectation <- mu + 2 ^ (1.0 / theta2) * sigma * gamma(2 / theta2) * (1 - 2 * theta1) * delta ^ (-1 / theta2) / gamma(1 / theta2)
    variance <- 2 ^ (2.0 / theta2) * sigma ^ 2 * (gamma(3 / theta2) * gamma(1 / theta2) * (1 - 3 * theta1 + 3 * theta1 ^ 2) - (gamma(2 / theta2)) ^ 2 * (1 - 2 * theta1) ^ 2) * delta ^ (-2 / theta2) / (gamma(1 / theta2)) ^ 2
    return(list(expectation = expectation, variance = variance))
}


  

