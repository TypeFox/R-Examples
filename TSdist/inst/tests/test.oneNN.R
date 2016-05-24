test_that("The that KMedoids function works properly", {
  
  # We create a training database with two series: 

  x <- cumsum(rnorm(50))
  y <- sin(seq(0, pi, length.out=50))
  
  train <- rbind(x, y)
  trainclasses <- c(1, 2)
  
  # We create a testing set with the same series, each replicated 5 times.
  
  test <- rbind(x, x, x, x, x, y, y, y, y, y)
  testclasses <- c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2)
  
  classification <- OneNN(train, trainclasses, test, testclasses, "euclidean")
  
  # Such easy database is always classified perfectly error=0. 
  expect_equal(classification$error, 0)
  
  # LCSS is a similarity, not a distance so, special treatment is necessary
  classification <- OneNN(train, trainclasses, test, testclasses, 
                          "lcss", epsilon=1)
  
  
  # If we don't define a ground truth classification for the testing set, 
  # we only get the class values
  
  classification <-OneNN(train, trainclasses, test, "euclidean")
  expect_equal(classification, testclasses)


  # If there are ties, the classes are selected randomly. The last series is just is at the same distance of x and y.
  test <- rbind(x, x, x, x, x, y, y, y, y, y, (x + y) /2)

  set.seed(123)
  classification <- OneNN(train, trainclasses, test, "euclidean")
  testclasses <- c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 1)
  expect_equal(classification, testclasses)

  set.seed(233)
  classification <- OneNN(train, trainclasses, test, "euclidean")
  testclasses <- c(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2)
  expect_equal(classification, testclasses)
  
  # If all the series in the training set are equal, special treatment 
  # is necessary
  train <- rbind(x, x)
  z <- y
  test <- rbind(z, z, z, z, z, z, z, z, z, z)
  classification <- OneNN(train, trainclasses, test, "euclidean")
  
})