context("phy.build")

test_that("congeneric.merge", {
  tree <- read.tree(text="((a_a:1,b_b:1):1, c_c:2):1;")
  tree <- congeneric.merge(tree, c("a_nother", "b_sharp"))
  expect_that(is.ultrametric(tree), is_true())
  expect_that(length(tree$tip.label), equals(5))
  tree <- read.tree(text="((a_a:1,b_b:1):1, c_c:2):1;")
  tree <- congeneric.merge(tree, c("a_nother", "a_gain", "b_sharp"))
  tree <- congeneric.merge(tree, c("a_nother", "b_sharp"))
  expect_that(is.ultrametric(tree), is_true())
  expect_that(length(tree$tip.label), equals(6))
})
