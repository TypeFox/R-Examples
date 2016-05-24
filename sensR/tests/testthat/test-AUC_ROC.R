context("AUC and ROC functions")

test_that("AUC works with negative d-primes", {
    ## Negative d-prime values allowed in AUC:
    res <- unname(unlist(lapply(c(-10, -1:1, 10), AUC)))
    res <- round(res, 7)
    ## dput(res)
    expect_equal(res, c(0, 0.2397501, 0.5, 0.7602499, 1))

    ## Also test confidence intervals:
    x <- lapply(-1:1, AUC, se.d=.2)
    ## dput(x)
    xRes <- list(structure(list(value = 0.239750061093477, lower = 0.162487075575793,
    upper = 0.333624729895317, CI.alpha = 0.05), .Names = c("value",
"lower", "upper", "CI.alpha"), class = "AUC"), structure(list(
    value = 0.5, lower = 0.390820654309087, upper = 0.609179345690913,
    CI.alpha = 0.05), .Names = c("value", "lower", "upper", "CI.alpha"
), class = "AUC"), structure(list(value = 0.760249938906523,
    lower = 0.666375270104683, upper = 0.837512924424207, CI.alpha = 0.05), .Names = c("value",
"lower", "upper", "CI.alpha"), class = "AUC"))

    expect_equal(x, xRes, tolerance=1e-7)
})

test_that("AUC works with scale argument", {
    dprime <- 5
    scale <- 2
    expect_equal(AUC(dprime, scale=scale)$value,
                 pnorm(dprime / sqrt(1 + scale^2)))
})


