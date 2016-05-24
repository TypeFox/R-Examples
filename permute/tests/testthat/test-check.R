library("testthat")
library("permute")

context("Testing check()")

## test that check will return all perms including the observed
test_that("check returns observed ordering in set of all permutations *if* asked to", {
    ## simple permutation
    h <- how(observed = TRUE)
    pp <- check(4, control = h)
    expect_that(nrow(pp$control$all.perms), equals(factorial(4)))

    ## time series
    h <- how(within = Within(type = "series"), observed = TRUE)
    n <- 10
    pp <- check(n, control = h)
    expect_that(nrow(pp$control$all.perms), equals(n))

    ## time series specified as a vector
    h <- how(within = Within(type = "series"), observed = TRUE)
    n <- 10
    vec <- seq_len(n)
    pp <- check(vec, control = h)
    expect_that(nrow(pp$control$all.perms), equals(n))
})

## test that check throws an error when within and plot permutation
## type is "none"
test_that("check detects if nothing to permute", {
    pl <- gl(4, 3)
    n <- 12
    h <- how(within = Within(type = "none"),
             plots = Plots(strata = pl, type = "none"))
    expect_error(check(seq_len(n), control = h),
                 regexp = "Permutation 'type' is \"none\" for both 'plots' & 'within'.\nNothing to permute.")
})
