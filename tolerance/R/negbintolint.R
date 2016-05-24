negbintol.int <- function (x, n, m = NULL, alpha = 0.05, P = 0.99, side = 1, method = c("LS", 
    "WU", "CB", "CS", "SC", "LR", "SP", "CC")) 
{
    if (side != 1 && side != 2) {
        stop(paste("Must specify a one-sided or two-sided procedure!", 
            "\n"))
    }
    if (side == 2) {
        alpha <- alpha/2
        P <- (P + 1)/2
    }
    if (is.null(m)) {
        m <- n
    }
    method <- match.arg(method)
    nu.hat <- n/(n + x)
    nu.tilde <- (n - 1)/(n + x - 1)
    z.a <- qnorm(1 - alpha)
    se.nu.hat <- sqrt((nu.hat^2 * (1 - nu.hat))/n)
    if (method == "LS") {
        lower.nu <- max(1e-07, nu.hat - z.a * se.nu.hat, na.rm = TRUE)
        upper.nu <- min(nu.hat + z.a * se.nu.hat, 1, na.rm = TRUE)
    }
    else if (method == "WU") {
        if ((n + x) <= 2) {
            stop(paste("Bounds are not defined for this option!", 
                "\n"))
        }
        TEMP <- nu.tilde + c(-1, 1) * z.a * sqrt(nu.tilde * (1 - 
            nu.tilde)/(n + x - 2))
        lower.nu <- max(1e-07, TEMP[1], na.rm = TRUE)
        upper.nu <- min(TEMP[2], 1, na.rm = TRUE)
    }
    else if (method == "CB") {
        TEMP <- c(n/(n + (x + 1) * qf(1 - alpha, 2 * (x + 1), 
            2 * n)), n * qf(1 - alpha, 2 * n, 2 * x)/(n * qf(1 - 
            alpha, 2 * n, 2 * x) + x))
        lower.nu <- max(1e-07, TEMP[1], na.rm = TRUE)
        upper.nu <- min(TEMP[2], 1, na.rm = TRUE)
    }
    else if (method == "CS") {
        TEMP <- qchisq(c(alpha, 1 - alpha), 2 * n)/(2 * (n + 
            x))
        lower.nu <- max(1e-07, TEMP[1], na.rm = TRUE)
        upper.nu <- min(TEMP[2], 1, na.rm = TRUE)
    }
    else if (method == "SC") {
        TEMP <- ((2 * (n + x) * n - n * z.a^2) + c(-1, 1) * sqrt(n^2 * 
            z.a^4 - 4 * (n + x) * n^2 * z.a^2 + 4 * (n + x)^2 * 
            n * z.a^2))/(2 * (n + x)^2)
        lower.nu <- max(1e-07, TEMP[1], na.rm = TRUE)
        upper.nu <- min(TEMP[2], 1, na.rm = TRUE)
    }
    else if (method == "LR") {
        fun.l = function(p, x, n, nu.hat, alpha) 2 * (log(dnbinom(x, 
            n, nu.hat)) - log(dnbinom(x, n, p))) - qchisq(1 - 
            alpha, 1)
        lower.nu <- suppressWarnings(try(uniroot(fun.l, c(1e-12, 
            nu.hat), x = x, n = n, nu.hat = nu.hat, alpha = alpha, 
            tol = 1e-20)$root, silent = TRUE))
        if (class(lower.nu) == "try-error") 
            lower.nu <- 1e-07
        fun.u = function(p, x, n, nu.hat, alpha) 2 * (log(dnbinom(x, 
            n, nu.hat)) - log(dnbinom(x, n, p))) - qchisq(alpha, 
            1)
        upper.nu <- suppressWarnings(try(uniroot(fun.u, c(nu.hat, 
            1), x = x, n = n, nu.hat = nu.hat, alpha = alpha, 
            tol = 1e-20)$root, silent = TRUE))
        if (class(upper.nu) == "try-error") 
            upper.nu <- 1
    }
    else if (method == "SP") {
        fun.sp <- function(p, x, n, K.a) {
            theta = log(x/((n + x) * (1 - p)))
            KY = n * (log(p) - log(1 - (1 - p) * exp(theta)))
            delta.1 = sign(theta) * sqrt(2 * abs(theta * x - 
                KY))
            delta.2 = theta * sqrt(x * (1 + x/n))
            FO = pnorm(delta.1) - dnorm(delta.1) * ((1/delta.2) - 
                (1/delta.1)) - K.a
            FO
        }
        lower.nu <- try(uniroot(fun.sp, c(1e-07, 0.9999999), 
            x = x, n = n, K.a = alpha, tol = 1e-20)$root, silent = TRUE)
        if (class(lower.nu) != "try-error") 
            if (abs(lower.nu - nu.hat) < 1e-05) 
                lower.nu <- try(uniroot(fun.sp, c(1e-07, 0.9999999), 
                  x = x, n = n, K.a = alpha, tol = 0.009)$root, 
                  silent = TRUE)
        if (class(lower.nu) == "try-error") 
            lower.nu <- 1e-07
        lower.nu <- max(1e-07, lower.nu, na.rm = TRUE)
        upper.nu <- try(uniroot(fun.sp, c(1e-12, 0.9999999), 
            x = x, n = n, K.a = 1 - alpha, tol = 1e-20)$root, 
            silent = TRUE)
        if (class(upper.nu) != "try-error") 
            if (abs(upper.nu - nu.hat) < 1e-05) 
                upper.nu <- try(uniroot(fun.sp, c(1e-12, 0.9999999), 
                  x = x, n = n, K.a = 1 - alpha, tol = 0.009)$root, 
                  silent = TRUE)
        if (class(upper.nu) == "try-error") 
            upper.nu <- 1
        upper.nu <- min(1, upper.nu, na.rm = TRUE)
        if (x == 0) {
            lower.nu <- 1
            upper.nu <- 1
        }
    }
    else if (method == "CC") {
        lower.nu <- max(1e-07, nu.hat - z.a * sqrt((nu.hat^2 * 
            (1 - nu.hat))/n) - 0.5/(n + x), na.rm = TRUE)
        upper.nu <- min(nu.hat + z.a * sqrt((nu.hat^2 * (1 - 
            nu.hat))/n) + 0.5/(n + x), 1, na.rm = TRUE)
    }
    lower <- qnbinom(1 - P, size = m, prob = upper.nu)
    if (lower.nu <= 1e-07) {
        upper <- Inf
    }
    else upper <- qnbinom(P, size = m, prob = lower.nu)
    if (side == 2) {
        alpha <- 2 * alpha
        P <- (2 * P) - 1
    }
    temp <- data.frame(cbind(alpha, P, nu.hat, lower, upper))
    if (side == 2) {
        colnames(temp) <- c("alpha", "P", "pi.hat", "2-sided.lower", 
            "2-sided.upper")
    }
    else {
        colnames(temp) <- c("alpha", "P", "pi.hat", "1-sided.lower", 
            "1-sided.upper")
    }
    temp
}
