context("cor_vector")

test_that("cor_vector only accepts numeric values", {
    m <- matrix(c(1L, 2L, 3L, 5L), ncol = 2)
    expect_error(cor_vector(m), "must be numeric")

    m <- matrix(c(1, 2, 3, 5), ncol = 2)
    expect_error(cor_vector(m), NA)
})

test_that("check case when all cors = 1", {
    m <- matrix(seq(1, by = 2, length.out = 12), ncol = 3)
    expect_equal(cor_vector(m), rep(1, 6))
})

test_that("cor_vector organizes values as expected", {
    m <- matrix(runif(12), ncol = 3) # random input here OK or no?
    expect_equal(cor_vector(m)[4:6], cor(m)[3,])
})

test_that("cor_vector has at least ncol ones", {
    m <- matrix(runif(80), ncol = 10) # random input here OK or no?
    expect_true(sum(cor_vector(m) == 1) >= 10)
})
