context("Perf")

test_that("Dataset with dissimilar categories", {
    set.seed(1)
    a <- cbind(1,rnorm(10,1,1))
    b <- cbind(5,rnorm(10,5,1))
    trainset <- as.data.frame(rbind(a,b))
    pnn <- learn(set=trainset)
    pnn$sigma <- 0.3
    pnn <- perf(pnn)
    expect_that(pnn$success_rate, equals(.9))
})

test_that("Dataset with letter categories", {
    set.seed(1)
    n <- 10
    m <- 1
    s <- 1
    d <- 4
    c <- rep(c("A","B"), each=n)
    x <- c( rnorm(n,m,s), rnorm(n,m,s) + d )
    trainset <- data.frame(c,x)
    pnn <- learn(trainset)
    pnn$sigma <- 0.3
    pnn <- perf(pnn)
    expect_that(pnn$success_rate, equals(.9))
})

test_that("Dataset with very similar categories", {
    set.seed(1)
    a <- cbind(1,rnorm(10,1,1))
    b <- cbind(2,rnorm(10,2,1))
    trainset <- as.data.frame(rbind(a,b))
    pnn <- learn(set=trainset)
    pnn$sigma <- 0.3
    pnn <- perf(pnn)
    expect_that(pnn$success_rate, equals(0.85))
    expect_that(pnn$bic, equals(-33.9, tolerance=0.1))
})
