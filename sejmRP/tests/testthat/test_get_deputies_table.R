test_that("result of function", {
  expect_equal(class(get_deputies_table(host = "services.mini.pw.edu.pl")), "data.frame")
})

test_that("columns of table", {
  expect_equal(ncol(get_deputies_table(host = "services.mini.pw.edu.pl")), 3)
})


test_that("rows of table", {
  expect_more_than(nrow(get_deputies_table(host = "services.mini.pw.edu.pl")), 0)
})