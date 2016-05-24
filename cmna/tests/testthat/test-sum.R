library("testthat")
context("sum")

l <- 1:10^6

n <- sample(l, 1)
bound <- sample(l, 2)
bound.u <- max(bound) - 10^6 / 2
bound.l <- min(bound) - 10^6 / 2
x <- runif(n, bound.l, bound.u)
expect_equal(naivesum(x), sum(x))
expect_equal(kahansum(x), sum(x))
expect_equal(pwisesum(x), sum(x))

n <- sample(l, 1)
bound <- sample(l, 2)
bound.u <- max(bound) - 10^6 / 2
bound.l <- min(bound) - 10^6 / 2
x <- runif(n, bound.l, bound.u)
expect_equal(naivesum(x), sum(x))
expect_equal(kahansum(x), sum(x))
expect_equal(pwisesum(x), sum(x))

n <- sample(l, 1)
bound <- sample(l, 2)
bound.u <- max(bound) - 10^6 / 2
bound.l <- min(bound) - 10^6 / 2
x <- runif(n, bound.l, bound.u)
expect_equal(naivesum(x), sum(x))
expect_equal(kahansum(x), sum(x))
expect_equal(pwisesum(x), sum(x))

n <- sample(l, 1)
bound <- sample(l, 2)
bound.u <- max(bound) - 10^6 / 2
bound.l <- min(bound) - 10^6 / 2
x <- runif(n, bound.l, bound.u)
expect_equal(naivesum(x), sum(x))
expect_equal(kahansum(x), sum(x))
expect_equal(pwisesum(x), sum(x))
