context("Mutability")

test_that("changing one affects copies", {
  a <- mutaframe(a = 1:10)
  b <- a
  
  a[1,1] <- 2
  expect_that(b[1, 1], equals(2))

  b[1,1] <- 3
  expect_that(a[1, 1], equals(3))
})

test_that("changing subset affects copies", {
  a <- mutaframe(a = 1:10)
  b <- a[1:5, , drop = FALSE]
  
  a[1,1] <- 2
  expect_that(b[1, 1], equals(2))
  expect_that(b[2:5, 1], equals(2:5))

  b[1,1] <- 3
  expect_that(a[1, 1], equals(3))
  expect_that(a[2:10, 1], equals(2:10))
})
