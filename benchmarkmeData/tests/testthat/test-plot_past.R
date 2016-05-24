test_that("Test plot_past", {
  skip_on_cran()
  expect_true(is.data.frame(plot_past()))
  expect_true(is.data.frame(plot_past(test=c("prog", "matrix_fun"))))
  expect_true(is.data.frame(plot_past(byte_optimize = FALSE)))
  expect_true(is.data.frame(plot_past(byte_optimize = TRUE)))
}
)
