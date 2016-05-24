
test_that("keeponly() returns the correct values", {
  vec1 <- c(0, 0, 0, 4, 4, 4, 100, 100, NA)
  vec2 <- c(0, 0.00000000001, 0, 4, 4, 4, 100, 100, 100)
  vec3 <- c(0, 0, 0, 4, 4, 4, 100, 100, 100.00000000001)
  vec4 <- c(0, 0, 0, 4, 4, 4, 100, 100, 100)
  vec5 <- c(0, 0, 4, 4, 4, 100, 100)
  vec6 <- c(0, 4, 4, 4, 100)
  vec7 <- c(4, 4, 4)

  expect_that(keeponly(vec1), throws_error())
  expect_that(keeponly(vec2), is_equivalent_to(c(FALSE, rep(TRUE, 7), FALSE)))
  expect_that(keeponly(vec3), is_equivalent_to(c(FALSE, rep(TRUE, 7), FALSE)))
  expect_that(keeponly(vec4), is_equivalent_to(c(FALSE, rep(TRUE, 7), FALSE)))
  expect_that(keeponly(vec5), is_equivalent_to(rep(TRUE, 7)))
  expect_that(keeponly(vec6), is_equivalent_to(rep(TRUE, 5)))
  expect_that(keeponly(vec7), is_equivalent_to(rep(TRUE, 3)))
})
