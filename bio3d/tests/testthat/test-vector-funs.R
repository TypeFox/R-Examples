context("Testing vector functions")


test_that("Vector functions", {

  ## vector normalization
  x <- 1:3
  x.norm <- c(0.2672612, 0.5345225, 0.8017837)
  expect_that(normalize.vector(x), equals(x.norm, tolerance  = 1e-6))

  y <- matrix(1:9, ncol = 3, nrow = 3)
  y.norm <- matrix(c(x.norm,
                     0.4558423, 0.5698029, 0.6837635,
                     0.5025707, 0.5743665, 0.6461623), ncol=3, byrow=F)

  expect_that(normalize.vector(y), equals(y.norm, tolerance  = 1e-6))


  
  ## Inner product
  x <- 1:3
  y <- diag(x)
  z <- matrix(1:9, ncol = 3, nrow = 3)

  xy <- c(1, 4, 9)
  yz <- c(1, 10, 27)
  
  expect_that(inner.prod(x,y), equals(xy))
  expect_that(inner.prod(y,z), equals(yz))
    
}
          )
