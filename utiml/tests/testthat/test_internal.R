context("R Internal tests")

test_that("Normalize", {
  expect_equal(utiml_normalize(1:3), c(0.0, 0.5, 1.0))
  expect_equal(utiml_normalize(1:5), c(0.0, 0.25, 0.5, 0.75, 1.0))
  expect_equal(utiml_normalize(c(1,2,3,4,5), 10, 0), 1:5/10)
})

test_that("Ifelse", {
  c1 <- rep(1, 10)
  c2 <- 1:10
  expect_equal(utiml_ifelse(T, c1, c2), c1)
  expect_equal(utiml_ifelse(F, c1, c2), c2)
  expect_equal(utiml_ifelse(T, c1, NA), c1)
  expect_equal(utiml_ifelse(F, c1, NA), NA)
  expect_equal(utiml_ifelse(T, NA, c2), NA)
  expect_equal(utiml_ifelse(F, NA, c2), c2)
  expect_null(utiml_ifelse(NA, c1, c2))
})

test_that("New data", {
  test <- toyml$dataset[,toyml$attributesIndexes]
  expect_equal(utiml_newdata(test), test)
  expect_equal(utiml_newdata(toyml), test)
})

test_that("Rename", {
  expect_equal(utiml_rename(c("a", "b", "c")), c(a="a", b="b", c="c"))
  expect_equal(utiml_rename(c(1, 2, 3), c("a", "b", "c")), c(a=1, b=2, c=3))
})

test_that("Equals sets", {
  expect_true(utiml_is_equal_sets(c(1, 2, 3), c(3, 2, 1)))
  expect_true(utiml_is_equal_sets(c(3, 2, 1), c(1, 2, 3)))
  expect_true(utiml_is_equal_sets(c(), c()))

  expect_false(utiml_is_equal_sets(c(1, 2, 3), c(1, 2, 3, 4)))
  expect_false(utiml_is_equal_sets(c(1, 2, 3, 4), c(1, 2, 3)))
  expect_false(utiml_is_equal_sets(c(1, 2, 3), c()))
})
