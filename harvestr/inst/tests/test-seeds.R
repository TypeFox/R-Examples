library(testthat)
context("seeds")

test_that("get.seed retrieves properly", {
    set.seed(123, "L'Ecuyer-CMRG", "Inversion")
    seed  <- get.seed()
    eseed <- c(407L
        , 1806547166L
        , -983674937L 
        ,  643431772L 
        , 1162448557L 
        , -959247990L 
        , -133913213L)
    expect_identical(seed, eseed)
})
test_that("replace.seed", {
    eseed <- c(407L
        , 1806547166L
        , -983674937L 
        ,  643431772L 
        , 1162448557L 
        , -959247990L 
        , -133913213L)
    set.seed(123, "Mersenne-Twister", "Box-Muller")
    k <- RNGkind()
    replace.seed(eseed)
    l <- RNGkind()
    
    expect_identical(k, c("Mersenne-Twister", "Box-Muller"))
    expect_identical(l, c("L'Ecuyer-CMRG", "Inversion"))
})
test_that("withseed, is replicable", {
    seeds <- gather(4)
    seed1 <- seeds[[1]]
    seed2 <- seeds[[2]]
    noattr(a <- withseed(seed1, rnorm(10)))
    noattr(b <- withseed(seed1, rnorm(10)))
    noattr(c <- withseed(seed2, rnorm(10)))
    d <- rnorm(10)
    expect_identical(a, b)
    expect_false(all(a == c))  # FALSE
    expect_false(identical(a, c))  # FALSE
    expect_false(all(a == d))  # FALSE
})
test_that("withseed replaces previous seed.", { 
    seed <- gather(1)[[1]]
    set.seed(123, "Mersenne-Twister", "Box-Muller")
    l <- get.seed()
    a <- withseed(seed, get.seed())
    k <- get.seed()
    expect_identical(l,k)
    expect_equivalent(a, attr(a, 'starting.seed'))
    expect_equivalent(a, attr(a, 'ending.seed'))
    expect_false(identical(l, noattr(a)))
    expect_false(isTRUE(all.equal(a, l)))
    expect_identical(RNGkind(), c("Mersenne-Twister", "Box-Muller"))
})
test_that('seeds set properly - withseed', {
    seeds <- gather(4)
    seed  <- noattr(seeds[[1]])
    
    set.seed(123, "Mersenne-Twister")
    cseed <- .Random.seed

    a <- withseed(seed, .Random.seed)

    s1 <- noattr(attr(a, 'starting.seed'))
    s2 <- noattr(a)
    s3 <- noattr(attr(a, 'ending.seed'))

    expect_identical(seed, s1)
    expect_identical(seed, s2)
    expect_identical(seed, s3)
    expect_false(isTRUE(all.equal(seed, cseed)))
    
    expect_identical(cseed, .Random.seed)
})
test_that('withseed handles passed expressions', {
    expr <- quote(runif(10))
    seed <- gather(1)[[1]]
    a <- withseed(seed, expr)
    b <- withseed(seed, runif(10))
    
    expect_that(a, equals(b))
    expect_that(length(a), equals(10))
    expect_that(all(0<=a & a<=1), is_true())
})
