

do.search <- function(phi, control, step0) {
    phi0 <- phi(0)
    ls <- linesearch(as.numeric(phi0), attr(phi0, "gradient"), step0, control)
    iter <- 0

    repeat {
        iter <- iter + 1
        phi1 <- phi(ls$step)
        ls <- update(ls, as.numeric(phi1), attr(phi1, "gradient"))

        if (ls$converged)
            break
    }

    list(iter = iter, step = ls$step, value = as.numeric(phi1),
         deriv = attr(phi1, "gradient"))
}


trunc <- function(x) {
    if (x == 0)
        return(0)
    if (x < 0)
        return(-trunc(-x))

    e <- floor(log10(x))

    if (e <= 1) {
        round(10^(-e+1) * x) / 10^(-e+1)
    } else {
        round(x / 10^(e-1)) * 10^(e-1)
    }
}


context("linesearch")

phi1 <- function(alpha, beta = 2) {
    res <- - alpha / (alpha^2 + beta)
    attr(res, "gradient") <- (alpha^2 - beta)/(alpha^2 + beta)^2
    res
}

control1 <- linesearch.control(value.tol = 0.001, deriv.tol = 0.1)

results1 <- rbind(c(1e-3, 6,  1.4, -9.2e-3),
                  c(1e-1, 3,  1.4,  4.4e-3), # reported 4.7e-3
                  c(1e+1, 1, 10  ,  9.4e-3),
                  c(1e+3, 4, 37  ,  7.3e-4))

test_that("Reproduce More and Thuente Table I", {
    for (i in seq_len(nrow(results1))) {
        actual <- do.search(phi1, control1, results1[i,1])
        expect_that(actual$iter, equals(results1[i,2]))
        expect_that(trunc(actual$step), equals(results1[i,3]))
        expect_that(trunc(actual$deriv), equals(results1[i,4]))
    }
})


phi2 <- function(alpha, beta = 0.004) {
    res <- (alpha + beta)^5 - 2 * (alpha + beta)^4
    attr(res, "gradient") <- 5 * (alpha + beta)^4  - 8 * (alpha + beta)^3
    res
}

control2 <- linesearch.control(value.tol = 0.1, deriv.tol = 0.1)

results2 <- rbind(c(1e-3, 12, 1.6,  3.8e-9),  # reported as 7.1e-9
                  c(1e-1,  8, 1.6,    1e-10), # reported as 10e-10; typo?
                  c(1e+1,  8, 1.6, -5.0e-9),
                  c(1e+3, 11, 1.6, -2.3e-8))

test_that("Reproduce More and Thuente Table II", {
    for (i in seq_len(nrow(results2))) {
        actual <- do.search(phi2, control2, results2[i,1])
        expect_that(actual$iter, equals(results2[i,2]))
        expect_that(trunc(actual$step), equals(results2[i,3]))
        expect_that(trunc(actual$deriv), equals(results2[i,4]))
    }
})


phi3 <- (function() {
    phi0 <- function(alpha, beta) {
        res <- ifelse(alpha <= 1 - beta, 1 - alpha,
                      ifelse(alpha >= 1 + beta, alpha - 1,
                             (alpha-1)^2 / (2 * beta) + beta/2))
        attr(res, "gradient") <-
            ifelse(alpha <= 1 - beta, -1,
                   ifelse(alpha >= 1 + beta, +1,
                          (alpha - 1) / beta))
        res
    }

    function(alpha, beta = 0.01, l = 39) {
        res0 <- phi0(alpha, beta)
        res <- (res0 + 2 * (1 - beta)/(l * pi) * sin(l * pi * alpha/2))
        attr(res, "gradient") <-
            (attr(res0, "gradient") + (1 - beta) * cos (l * pi / 2 * alpha))
        res
    }
})()

control3 <- linesearch.control(value.tol = 0.1, deriv.tol = 0.1)

results3 <- rbind(c(1e-3, 12, 1.0, -9.1e-5), # reported -5.1e-5
                  c(1e-1, 11, 1.0,  9.1e-7), # reported 12, -1.9e-4
                  c(1e+1,  9, 1.0, -5.3e-4), # reported 10, -2.0e-6
                  c(1e+3, 11, 1.0,  9.8e-6)) # reported 13, -1.6e-5

test_that("Reproduce More and Thuente Table III", {
    for (i in seq_len(nrow(results3))) {
        actual <- do.search(phi3, control3, results3[i,1])
        expect_that(actual$iter, equals(results3[i,2]))
        expect_that(trunc(actual$step), equals(results3[i,3]))
        expect_that(trunc(actual$deriv), equals(results3[i,4]))
    }
})


yanai <- function(beta1, beta2) {
    gamma <- function(beta) {
        sqrt(1 + beta^2) - beta
    }

    function(alpha) {
        res <- (gamma(beta1) * sqrt((1 - alpha)^2 + beta2^2)
                + gamma(beta2) * sqrt(alpha^2 + beta1^2))
        attr(res, "gradient") <-
            (gamma(beta1) * (alpha - 1) / sqrt((1 - alpha)^2 + beta2^2)
             + gamma(beta2) * alpha / sqrt(alpha^2 + beta1^2))
        res
    }
}


phi4 <- yanai(0.001, 0.001)
control4 <- linesearch.control(value.tol = 0.001, deriv.tol = 0.001)
results4 <- rbind(c(1e-3, 4, 0.085, -6.9e-5), # reported as 0.08; rounding difference?
                  c(1e-1, 1, 0.10,  -4.9e-5),
                  c(1e+1, 3, 0.34,  -3.2e-6), # reported 0.35, -2.9e-6
                  c(1e+3, 4, 0.83,   1.6e-5))

test_that("Reproduce More and Thuente Table IV", {
    for (i in seq_len(nrow(results4))) {
        actual <- do.search(phi4, control4, results4[i,1])
        expect_that(actual$iter, equals(results4[i,2]))
        expect_that(trunc(actual$step), equals(results4[i,3]))
        expect_that(trunc(actual$deriv), equals(results4[i,4]))
    }
})


phi5 <- yanai(0.01, 0.001)
control5 <- linesearch.control(value.tol = 0.001, deriv.tol = 0.001)
results5 <- rbind(c(1e-3, 6, 0.075,  1.9e-4),
                  c(1e-1, 3, 0.078,  7.4e-4),
                  c(1e+1, 7, 0.073, -2.5e-4), # reported -2.6e-4
                  c(1e+3, 8, 0.076,  4.4e-4)) # reported 4.5e-4

test_that("Reproduce More and Thuente Table V", {
    for (i in seq_len(nrow(results5))) {
        actual <- do.search(phi5, control5, results5[i,1])
        expect_that(actual$iter, equals(results5[i,2]))
        expect_that(trunc(actual$step), equals(results5[i,3]))
        expect_that(trunc(actual$deriv), equals(results5[i,4]))
    }
})


phi6 <- yanai(0.001, 0.01)
control6 <- linesearch.control(value.tol = 0.001, deriv.tol = 0.001)
results6 <- rbind(c(1e-3, 13, 0.93,  5.6e-4), # reported as 5.2e-4
                  c(1e-1, 11, 0.93,  2.3e-4), # reported as 8.4e-5
                  c(1e+1,  8, 0.92, -2.1e-4), # reported as -2.4e-4
                  c(1e+3, 10, 0.92, -3.7e-4)) # reported as 11 and -3.2e-4

test_that("Reproduce More and Thuente Table VI", {
    for (i in seq_len(nrow(results6))) {
        actual <- do.search(phi6, control6, results6[i,1])
        expect_that(actual$iter, equals(results6[i,2]))
        expect_that(trunc(actual$step), equals(results6[i,3]))
        expect_that(trunc(actual$deriv), equals(results6[i,4]))
    }
})

