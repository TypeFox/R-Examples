context("gen_lambdas")

test_that("gen_lambdas is the same as generate.lambdas", {
    lambda.max <- 10
    lambda.min <- 1
    expect_equal(gen_lambdas(lambda.max, lambda.min), generate.lambdas(lambda.max, lambda.min / lambda.max))
})

test_that("gen_lambdas produces values in the expected range", {
    max.val <- sqrt(10)
    min.val <- 0.1*max.val

    ### Linear scale
    expect_true(all(gen_lambdas(max.val, min.val, scale = "linear") >= 0))
    expect_true(all(gen_lambdas(max.val, min.val, scale = "linear") <= max.val))
    expect_equal(max(gen_lambdas(max.val, min.val, scale = "linear")), max.val)
    expect_equal(min(gen_lambdas(max.val, min.val, scale = "linear")), min.val)

    ### Log scale
    expect_true(all(gen_lambdas(max.val, min.val, scale = "log") >= 0))
    expect_true(all(gen_lambdas(max.val, min.val, scale = "log") <= max.val))
    expect_equal(max(gen_lambdas(max.val, min.val, scale = "log")), max.val)
    expect_equal(min(gen_lambdas(max.val, min.val, scale = "log")), min.val)
})

test_that("gen_lambdas produces the exact output expected on test cases", {
    vec <- gen_lambdas(10, 1, lambdas.length = 10, scale = "linear") # should be 10:1
    expect_equal(vec, 10:1)

    vec <- gen_lambdas(10, 1, lambdas.length = 3, scale = "log") # should be {10.000000, 3.162278, 1.000000}
    expect_true(sum(abs(vec - c(10.000000, 3.162278, 1.000000))) < 1e-6) # need to check error instead of equality since values are approximate
})
