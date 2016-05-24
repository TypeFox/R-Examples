test_that("Testing get_all_typeforms", {
  skip_on_cran()
  typeforms = get_all_typeforms()
  expect_equal(ncol(typeforms), 2)
}
)