test_that("Test datatable_past", {
  skip_on_cran()
  expect_true(is.list(get_datatable_past()))
  expect_true(is.list(get_datatable_past(test=c("prog", "matrix_fun"))))
  expect_true(is.list(get_datatable_past(byte_optimize = TRUE)))
  expect_true(is.list(get_datatable_past(byte_optimize = FALSE)))
  }
)
