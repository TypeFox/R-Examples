context("Smooth")

test_that("Set sigma", {
    data(norms)
    pnn <- learn(norms)
    pnn <- smooth(pnn, sigma=0.8)
    expect_that(pnn$sigma, equals(0.8))
})

test_that("Function to minimize", {
    set.seed(1)
    a <- cbind(1,rnorm(10,1,1))
    b <- cbind(2,rnorm(10,2,1))
    trainset <- as.data.frame(rbind(a,b))
    pnn <- learn(set=trainset)
    tominimize <- iToMinimize(pnn)
    expect_that(is.function(tominimize), is_true())
    expect_that(tominimize(0.3), equals(3))
})