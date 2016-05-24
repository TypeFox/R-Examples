context("firthglm.fit")

test_that("succeeds on dataset 1", {
    data <- readRDS("firth-binomial1.rds")
    model <- firthglm.fit(data$x, data$y, family=binomial())
    coef <- c(11.4, 12.7, -6.9, 0.7, 14.8, -2.4)

    expect_that(model$converged, is_true())
    expect_that(model$boundary, is_false())
    expect_that(as.numeric(round(model$coefficients, 1)), equals(coef))
})


test_that("succeeds on dataset 2", {
    data <- readRDS("firth-binomial2.rds")
    model <- firthglm.fit(data$x, data$y, family=binomial())
    coef <- c(-10.3, 1.6, -1.3, 0.5, 3.1, 0.8)

    expect_that(model$converged, is_true())
    expect_that(model$boundary, is_false())
    expect_that(as.numeric(round(model$coefficients, 1)), equals(coef))
})


test_that("succeeds on dataset 3", {
    data <- readRDS("firth-binomial3.rds")
    model <- firthglm.fit(data$x, data$y, family=binomial())
    coef <- c(12.1, 0.0, -2.9, NA, 0.1, -1.1)

    expect_that(model$converged, is_true())
    expect_that(model$boundary, is_false())
    expect_that(as.numeric(round(model$coefficients, 1)), equals(coef))
})


test_that("succeeds on dataset 4", {
    data <- readRDS("firth-binomial4.rds")
    model <- firthglm.fit(data$x, data$y, family=binomial())
    coef <- c(-100.8, 18.2, 2.3, NA, 4.0, 12.6)

    expect_that(model$converged, is_true())
    expect_that(model$boundary, is_false())
    expect_that(as.numeric(round(model$coefficients, 1)), equals(coef))
})


test_that("succeeds on dataset 5", {
    data <- readRDS("firth-binomial5.rds")
    model <- firthglm.fit(data$x, data$y, family=binomial())
    coef <- c(-6.0, 7.3, 21.7, 4.9, 7.6, 10.6, 7.0, 9.1, 5.4, -0.7, 4.9)

    expect_that(model$converged, is_true())
    expect_that(model$boundary, is_false())
    expect_that(as.numeric(round(model$coefficients, 1)), equals(coef))
})


test_that("succeeds on dataset 6", {
    data <- readRDS("firth-binomial6.rds")
    model <- firthglm.fit(data$x, data$y, family=binomial())
    coef <- c(-1.1, -0.2, -0.4, 0.1,  0.4, -0.1, -1.4,  0.7,  0.8, 0.9,
              -0.4,  0.6, -0.2, 0.2, -0.6,  0.2,  2.3, -1.7, -1.2)

    expect_that(model$converged, is_true())
    expect_that(model$boundary, is_false())
    expect_that(as.numeric(round(model$coefficients, 1)), equals(coef))
})

test_that("succeeds on dataset 7", {
    data <- readRDS("firth-binomial7.rds")
    model <- firthglm.fit(data$x, data$y, family=binomial())
    coef <- c(-7.1, 2.1, 5.7)

    expect_that(model$converged, is_true())
    expect_that(model$boundary, is_false())
    expect_that(as.numeric(round(model$coefficients, 1)), equals(coef))
})

test_that("succeeds on dataset 8", {
    data <- readRDS("firth-binomial8.rds")
    model <- firthglm.fit(data$x, data$y, family=binomial())
    coef <- c(-24.7, -12.4, 12.0, 2.2, 19.3, -16.6, -6.0, -8.1, 25.0, 22.5)

    expect_that(model$converged, is_true())
    expect_that(model$boundary, is_false())
    expect_that(as.numeric(round(model$coefficients, 1)), equals(coef))
})

test_that("succeeds on dataset 9", {
    data <- readRDS("firth-binomial9.rds")
    model <- firthglm.fit(data$x, data$y, family=binomial(), start=data$start)
    coef <- c( 0.0,  0.1,  0.3, -0.2, -0.2,  0.1, -0.3,  0.2,  0.2, -0.1,
               0.0,  0.0,  0.0,  0.1,  0.2,  0.0,  0.1, -0.2,  0.1,  0.2,
              -0.3,  0.0, -0.2, -0.3,  0.0, -0.1,  0.3, -0.1, -0.2, -0.3,
               0.1,  0.5, -0.2, -0.1,  0.4, -0.1,  0.0, -0.1,  0.2,  0.2,
               0.3, -0.1,  0.3,  0.3,  0.1,  0.1, -0.2,  0.0,  0.0,  0.5)

    expect_that(model$converged, is_true())
    expect_that(model$boundary, is_false())
    expect_that(as.numeric(round(model$coefficients, 1)), equals(coef))
})

test_that("succeeds on dataset 10", {
    data <- readRDS("firth-binomial10.rds")
    model <- firthglm.fit(data$x, data$y, family=binomial())
    coef <- c(-0.1,  0.1, 0.2,  0.2, -0.1,  0.3,  0.1,  0.0,  0.0,  0.0,
              -0.1,  0.1, 0.2, -0.1,  0.0,  0.2, -0.2, -0.2,  0.1,  0.1,
              -0.1, -0.3, 0.5,  0.0, -0.1,  0.2,  0.2,  0.0,  0.1,  0.3,
              -0.1,  0.0, 0.0, -0.3,  0.1, -0.1,  0.0, -0.2, -0.3, -0.1,
              -0.1,  0.0, 0.3,  0.2,  0.1, -0.1,  0.4,  0.0,  0.0,  0.0,
               0.1,  0.2, 0.2, -0.2,  0.0,  0.0,  0.1,   NA,   NA,   NA,
                NA,   NA,  NA,   NA,   NA,   NA,   NA,   NA,   NA,   NA,
                NA,   NA,  NA,   NA,   NA,   NA,   NA,   NA,   NA,   NA,
                NA,   NA,  NA,   NA,   NA,   NA,   NA,   NA,   NA,   NA,
                NA,   NA,  NA,   NA,   NA,   NA,   NA,   NA,   NA,   NA)

    expect_that(model$converged, is_true())
    expect_that(model$boundary, is_false())
    expect_that(as.numeric(round(model$coefficients, 1)), equals(coef))
})

test_that("succeeds on rank degenerate case", {
    x <- matrix(c(-1, 1, -1, 1, 1, 1, 1, -1, 1, -1), nrow=2)
    y <- c(1, 1)

    model <- firthglm.fit(x, y, family=binomial())
    expect_that(model$converged, is_true())
})
