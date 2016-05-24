context("checkListElementClass")

test_that("checkListElementClass", {
  checkListElementClass(list(1, 5), cl="numeric")
  expect_error(checkListElementClass(list(1, "a"), cl="numeric"), "numeric")
  
  xs = list("a", "b")
  checkListElementClass(xs, "character")
  expect_error(checkListElementClass(xs, "integer"), "character")
})


