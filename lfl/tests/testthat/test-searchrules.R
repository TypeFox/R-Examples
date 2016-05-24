test_that('searchrules', {
    n <- 100
    d <- data.frame(a=1:n, b=n:1, c=runif(n))
    d <- lcut3(d)

    rules <- searchrules(d)
    expect_true(is.farules(rules))
})
