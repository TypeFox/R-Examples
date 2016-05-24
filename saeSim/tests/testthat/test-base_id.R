context("base_id: ")
test_that("base_id", {
  expect_equal(
    base_id(3, c(1, 2, 3)), 
    data.frame(idD = c(1, 2, 2, 3, 3, 3), idU = c(1, 1, 2, 1, 2, 3))
  )
  expect_equal(base_id(3, 1), data.frame(idD = c(1, 2, 3)))
})

test_that("base_id_temporal", {
  
  dat <- base_id_temporal(2, 2, 2)
  
  comp <- data.frame(
    idD = rep(1:2, each = 4), 
    idU = c(1, 1, 2, 2, 1, 1, 2, 2), 
    idT = rep(1:2, 4)
  )
  
  testthat::expect_equal(dat, comp, check.attributes = FALSE)
  
  dat <- base_id_temporal(2, 1, 2)
  
  comp <- data.frame(
    idD = rep(1:2, each = 2), 
    idT = rep(1:2, 2)
  )
  
  testthat::expect_equal(dat, comp, check.attributes = FALSE)
  
  dat <- base_id_temporal(2, 1:2, 2)
  
  comp <- data.frame(
    idD = c(1, 1, 2, 2, 2, 2), 
    idU = c(1, 1, 1, 1, 2, 2),
    idT = c(1, 2, 1, 2, 1, 2)
  )
  
  testthat::expect_equal(dat, comp, check.attributes = FALSE)
  
})