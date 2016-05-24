test_that("The that KMedoids function works properly", {
  
  # We create a database with two series, each one replicated 5 times: 

  x <- cumsum(rnorm(50))
  y <- sin(seq(0, pi, length.out=50))
  
  data <- rbind(x, x, x, x, x, y, y, y, y, y)
  classes <- c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2)
  
  clustering <- KMedoids(data, k=2, classes, "euclidean")
  
  #Such easy database is always clustered perfectly F=1. 
  expect_equal(clustering$F, 1)
  
  #If we don't define a ground truth, we only get the clustering value
  clustering <- as.numeric(KMedoids(data, k=2, "euclidean"))
  expect_equal(clustering, classes)
  
  # LCSS is not a distance measure, but a similarity, so a different treatment
  # is necessary 
  clustering <- as.numeric(KMedoids(data, k=2, "lcss", epsilon=1))
  expect_equal(clustering, classes)


})
