context("Basic functionality")


test_that("character method functions as expected", {
  
  m1 <- mat("1.5, 2.0, pi; 4, 5, 6; 7, 8, 9")
  m2 <- mat("3/2 2.0 pi  ; 4 5 6;    7:9 ", sep = " ")
  m3 <- matrix(c(1.5, 2, pi, 4:9), nrow = 3, byrow = TRUE)
  expect_that(m1, is_identical_to(m2))
  expect_that(m1, equals(m3, check.attributes = FALSE))
  expect_that(m2, equals(m3, check.attributes = FALSE))
  
})


test_that("list method functions as expected", {
  
  z <- list(a = 1:10, b = 1:10, c = 1:10)
  m1 <- mat(z)
  m2 <- mat(z, rows = FALSE)
  m3 <- matrix(unlist(z), nrow = 3, byrow = TRUE)
  m4 <- do.call(rbind, z)
  m5 <- simplify2array(z)
  expect_that(rownames(m1), is_identical_to(names(z)))  # check names
  expect_that(colnames(m2), is_identical_to(names(z)))  # check names
  expect_that(m1, is_identical_to(t(m2)))
  expect_that(m1, equals(m4, check.attributes = FALSE))
  expect_that(m2, equals(m5, check.attributes = FALSE))
  expect_that(m1, equals(m3, check.attributes = FALSE))
  
})


test_that("dmat functions as expected", {
  
  z <- list(a = 1:10, b = 1:10, c = 1:10)
  d1 <- dmat(z)
  d2 <- dmat(z, rows = FALSE)
  m1 <- mat(z)
  m2 <- mat(z, rows = FALSE)
  expect_that(d1, is_identical_to(as.data.frame(m1)))
  expect_that(d2, is_identical_to(as.data.frame(m2)))
  expect_that(m1, is_identical_to(t(m2)))
  
})