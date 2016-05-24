test_that("package loads successfully", {
  data("treebase")
  expect_more_than(length(treebase), 1)
})
