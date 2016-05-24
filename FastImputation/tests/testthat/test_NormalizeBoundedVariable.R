context("NormalizeBoundedVariable")

test_that("NormalizeBoundedVariable catches bad input", {
  good.x.vector <- c(-3:3)
  good.x.matrix <- matrix(-2:3, nrow=2)
  good.x.df <- data.frame(x=-3:0, y=0:3)
  good.x.array <- array(c(-3:0, 0:3), dim=c(4,2))
  
  bad.x.list <- list(x=1:3, y=4:6)
  
  bad.constraints.order <- list(lower=1, upper=0)
  good.constraints.lower <- list(lower=-5)
  good.constraints.upper <- list(upper=5)
  good.constraints.both <- list(lower=-5, upper=5)  
  
  expect_error(BoundNormalizedVariable(x=good.x.vector, constraints=bad.constraints.order), 
    "'upper' must be greater than 'lower.'")
  expect_error(BoundNormalizedVariable(x=good.x.matrix, constraints=bad.constraints.order), 
    "'upper' must be greater than 'lower.'")
  expect_error(BoundNormalizedVariable(x=good.x.df, constraints=bad.constraints.order), 
    "'upper' must be greater than 'lower.'")
  expect_error(BoundNormalizedVariable(x=good.x.array, constraints=bad.constraints.order), 
    "'upper' must be greater than 'lower.'")
})

test_that("NormalizeBoundedVariable returns correct types", {
  good.x.vector <- c(-3:3)
  good.x.matrix <- matrix(-2:3, nrow=2)
  good.x.df <- data.frame(x=-3:0, y=0:3)
  good.x.array <- array(c(-3:0, 0:3), dim=c(4,2))
  good.constraints.lower <- list(lower=-5)
  good.constraints.upper <- list(upper=5)
  good.constraints.both <- list(lower=-5, upper=5)  
  
  expect_true( is.vector(BoundNormalizedVariable(x=good.x.vector, 
    constraints=good.constraints.lower)) )
  expect_true( is.vector(BoundNormalizedVariable(x=good.x.vector, 
    constraints=good.constraints.upper)) )
  expect_true( is.vector(BoundNormalizedVariable(x=good.x.vector, 
    constraints=good.constraints.both)) )

  expect_true( is.matrix(BoundNormalizedVariable(x=good.x.matrix, 
    constraints=good.constraints.lower)) )
  expect_true( is.matrix(BoundNormalizedVariable(x=good.x.matrix, 
    constraints=good.constraints.upper)) )
  expect_true( is.matrix(BoundNormalizedVariable(x=good.x.matrix, 
    constraints=good.constraints.both)) )

  expect_true( is.data.frame(BoundNormalizedVariable(x=good.x.df, 
    constraints=good.constraints.lower)) )
  expect_true( is.data.frame(BoundNormalizedVariable(x=good.x.df, 
    constraints=good.constraints.upper)) )
  expect_true( is.data.frame(BoundNormalizedVariable(x=good.x.df, 
    constraints=good.constraints.both)) )

  expect_true( is.array(BoundNormalizedVariable(x=good.x.array, 
    constraints=good.constraints.lower)) )
  expect_true( is.array(BoundNormalizedVariable(x=good.x.array, 
    constraints=good.constraints.upper)) )
  expect_true( is.array(BoundNormalizedVariable(x=good.x.array, 
    constraints=good.constraints.both)) )
})

test_that("NormalizeBoundedVariable returns correct values", {
  good.x <- 1
  good.constraints.lower <- list(lower=-5)
  good.constraints.upper <- list(upper=5)
  good.constraints.both <- list(lower=-5, upper=5)  
  
  expect_equal( BoundNormalizedVariable(x=good.x, constraints=good.constraints.lower), 
    exp(good.x) + good.constraints.lower$lower )
  expect_equal( BoundNormalizedVariable(x=good.x, constraints=good.constraints.upper),
    good.constraints.upper$upper - exp(good.x) )
  expect_equal( BoundNormalizedVariable(x=good.x, constraints=good.constraints.both), 
    pnorm(good.x) * (good.constraints.both$upper - good.constraints.both$lower) + 
      good.constraints.both$lower )
      
  #expect_true(FALSE)
})
