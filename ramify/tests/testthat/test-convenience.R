context("Convenience functions")


test_that("convenience functions work as expected", {
  
  # argmax/argmin
  a <- resize(0:5, 2, 3)
  expect_identical(argmax(a), c(3L, 3L))
  expect_identical(argmax(a, rows = FALSE), c(2L, 2L, 2L))
  expect_identical(argmin(a), c(1L, 1L))
  expect_identical(argmin(a, rows = FALSE), c(1L, 1L, 1L))
  
  # eye
  d <- eye(3, 3)
  expect_identical(d, diag(1L, 3, 3))
  expect_identical(d, inv(d))  # inverse of identity is identity
  
  # fill (ones, zeros, trues, falses), rand, randi, randn
  expect_equal(dim(fill(0L, 10)), c(10, 1))
  expect_equal(dim(randn(10)), c(10, 1))
  expect_equal(dim(rand(10)), c(10, 1))
  expect_equal(dim(randi(imax = 100, 10)), c(10, 1))
  expect_null(dim(fill(0L, 10, atleast_2d = FALSE)))
  expect_null(dim(randn(10, atleast_2d = FALSE)))
  expect_null(dim(rand(10, atleast_2d = FALSE)))
  expect_null(dim(randi(imax = 100, 10, atleast_2d = FALSE)))
  expect_that(dim(fill(0L, 2, 2, 2)), equals(c(2, 2, 2)))
  expect_that(dim(rand(2, 2, 2)), equals(c(2, 2, 2)))
  expect_that(dim(randi(imax = 100, 2, 2, 2)), equals(c(2, 2, 2)))
  expect_that(dim(randn(2, 2, 2)), equals(c(2, 2, 2)))
  
  # flatten
  m1 <- matrix(1:9, 3, 3, byrow = TRUE)
  m2 <- matrix(1:9, 3, 3, byrow = FALSE)
  expect_identical(flatten(m1), 1:9)
  expect_identical(flatten(m1), flatten(m2, across = "columns"))
  expect_identical(flatten(m2), as.integer(c(1, 4, 7, 2, 5, 8, 3, 6, 9)))
  
  # inv, tr
  expect_error(inv(fill(3L, 2, 2)))  # singular matrix
  expect_identical(inv(2 * eye(3, 3)), 0.5 * eye(3, 3))
  expect_identical(tr(5 * eye(3, 4)), 15)
  
  # hcat, vcat
  m1 <- randn(2, 3)
  m2 <- randn(2, 3)
  expect_identical(hcat(m1, m2), cbind(m1, m2))
  expect_identical(vcat(m1, m2), rbind(m1, m2))
  
  # linspace, logspace
  
  # meshgrid
  x <- linspace(0, 1, 3)
  y <- linspace(0, 1, 2)
  expect_identical(meshgrid(x, y), list(mat("0, 0.5, 1; 0, 0.5, 1"),
                                        mat("0, 0,   0; 1, 1,   1")))
  x <- y <- seq(-5, 5, by = 0.1)
  mg <- meshgrid(x, y)
  z1 <- sin(mg[[1]]^2 + mg[[2]]^2) / (mg[[1]]^2 + mg[[2]]^2)
  z2 <- outer(x, y, function(x, y) sin(x^2 + y^2) / (x^2 + y^2))
  expect_identical(z1, z2)
  
  # resize, size
  
  # tri, tril, triu, is.tril, is.triu

  # Resize a vector into an array
  x <- 1:8
  a <- resize(1:8, 2, 2, 2)
  expect_that(a, is_a("array"))
  expect_that(flatten(a), is_identical_to(x))
  

  # Meshgrid (identical)
  x <- y <- seq(-5, 5, by = 0.1)
  mg <- meshgrid(x, y)
  z1 <- sin(mg[[1]]^2 + mg[[2]]^2) / (mg[[1]]^2 + mg[[2]]^2)
  z2 <- outer(x, y, function(x, y) sin(x^2 + y^2) / (x^2 + y^2))
  expect_that(z1, is_identical_to(z2))
  
  # Triangular matrices
  m1 <- mat("1, 1, 1, 0, 0; 
             1, 1, 1, 1, 0; 
             1, 1, 1, 1, 1")
  m2 <- mat("0, 0, 0, 0, 0; 
             1, 0, 0, 0, 0; 
             1, 1, 0, 0, 0")
  expect_that(tri(3, 5, k = 2), equals(m1, check.attributes = FALSE))
  expect_that(tri(3, 5, k = -1), equals(m2, check.attributes = FALSE))                       
  expect_that(tri(3, 5, diag = FALSE), equals(m2, check.attributes = FALSE))
  
  m3 <- matrix(1:12, nrow = 4, ncol = 3, byrow = TRUE)
  m4 <- mat(" 0,  0,  0; 
              4,  0,  0; 
              7,  8,  0; 
             10, 11, 12")
  expect_that(tril(m3, k = -1), equals(m4, check.attributes = FALSE))
  
})
