test_that('fire antecedents on data vector', {
    x <- 1:10 / 10
    names(x) <- letters[1:10]
    rules <- list(c('a', 'c', 'e'),
                  c('b'),
                  c('d', 'a'),
                  c('c', 'a', 'b'))

    res <- fire(x, rules, goguen.tnorm)
    expect_equal(res, list(0.15, 1, 0.1, 0.02))

    res <- fire(x, rules, 'goguen')
    expect_equal(res, list(0.15, 1, 0.1, 0.02))
})


test_that('fire single antecedent on data matrix', {
    x <- matrix(1:20 / 20, nrow=2)
    colnames(x) <- letters[1:10]
    rules <- c('a', 'c', 'e')

    res <- fire(x, rules, goguen.tnorm)

    expect_equal(res, list(c(0.1125, 0.15)))
})


test_that('fire whole single rule on data matrix', {
    x <- matrix(1:20 / 20, nrow=2)
    colnames(x) <- letters[1:10]
    rules <- c('a', 'c', 'e')

    res <- fire(x, rules, goguen.tnorm, onlyAnte=FALSE)

    expect_equal(res, list(c(0.005625, 0.015)))
})


test_that('fire antecedents on data matrix', {
    x <- matrix(1:20 / 20, nrow=2)
    colnames(x) <- letters[1:10]
    rules <- list(c('a', 'c', 'e'),
                  c('b'),
                  c('d', 'a'),
                  c('c', 'a', 'b'))

    res <- fire(x, rules, goguen.tnorm)

    expect_equal(res, list(c(0.1125, 0.15),
                           c(1, 1),
                           c(0.05, 0.1),
                           c(0.0075, 0.02)))
})


test_that('fire on farules', {
    x <- matrix(1:20 / 20, nrow=2)
    colnames(x) <- letters[1:10]
    rules <- list(c('a', 'c', 'e'),
                  c('b'),
                  c('d', 'a'),
                  c('c', 'a', 'b'))
    farules <- list(rules=rules, statistics=matrix(0, ncol=1, nrow=1))
    class(farules) <- c('farules', 'list')
    expect_true(is.farules(farules))

    res <- fire(x, farules, goguen.tnorm)

    expect_equal(res, list(c(0.1125, 0.15),
                           c(1, 1),
                           c(0.05, 0.1),
                           c(0.0075, 0.02)))
})
