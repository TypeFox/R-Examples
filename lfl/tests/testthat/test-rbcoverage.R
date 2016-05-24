test_that('rbcoverage on data vector', {
    x <- 1:10 / 10
    names(x) <- letters[1:10]
    rules <- list(c('a', 'c', 'e'),
                  c('b'),
                  c('d', 'a'),
                  c('c', 'a', 'b'))

    expect_equal(rbcoverage(x, rules, "goguen", TRUE),
                 1)
})


test_that('rbcoverage on data matrix', {
    x <- matrix(1:20 / 20, nrow=2)
    colnames(x) <- letters[1:10]

    rules <- list(c('a', 'c', 'e'),
                  c('b'),
                  c('d', 'a'),
                  c('c', 'a', 'b'))
    expect_equal(rbcoverage(x, rules, "goguen", TRUE),
                 1)

    rules <- list(c('a', 'c', 'e'),
                  c('d', 'a'),
                  c('c', 'a', 'b'))
    expect_equal(rbcoverage(x, rules, "goguen", TRUE),
                 0.13125)

    rules <- list(c('d', 'a'),
                  c('c', 'a', 'b'))
    expect_equal(rbcoverage(x, rules, "goguen", TRUE),
                 0.075)

    rules <- list(c('d', 'a'))
    expect_equal(rbcoverage(x, rules, "goguen", TRUE),
                 0.075)
})

