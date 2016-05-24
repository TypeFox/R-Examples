context("Testing p_m1\n")

test_that("p_m1 for N = 1 equals P(U < 1/gamma)", {
  result <- p_m1(0.1, beta = c(2, 3), "normal")
  expect_equal(result, pnorm(-1 / 0.1))
})
