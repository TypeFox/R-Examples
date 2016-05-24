test_that("Testing get_results_typeforms", {
  skip_on_cran()
  uid = "o99fIa"
  res = get_results(uid)
  expect_true(is.data.frame(res))
}
)