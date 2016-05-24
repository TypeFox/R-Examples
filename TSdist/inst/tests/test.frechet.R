context("Frechet")

test_that("The Frechet distance is calculated correctly", {
  
  # We use an example in the documentation of the longitudinalData package to test that the wrapper works correctly: 
  
  Px <- 1:20
  Py <- dnorm(1:20, 12, 2)
  Qx <- 1:20
  Qy <- dnorm(1:20, 8, 2)
  
  d<-distFrechet(Px,Py,Qx,Qy)

  expect_equal(FrechetDistance(x=Py, y=Qy, tx=Px, ty=Qx), d)
  
  #If the original function outputs an error, the wrapper will return NA
  expect_equal(FrechetDistance(x=Py, y=Qy, tx=Px, ty=Qx, timeScale="a"), NA)

  
})
