context("TQuest")

test_that("The TQuest function is calculated correctly", {
  
  # FIrst example
  x <- c(5, 5, 7, 7, 7, 5, 5)
  y <- c(2, 2, 8, 8, 8, 2, 2)
  tx <- c(1, 2, 3, 4, 5, 6, 7)
  ty <- c(1, 2, 3, 4, 5, 6, 7)

  
  # If we set tau=6, both these series trespass this threshold in time 
  # intervals {3,4,5}. So the distance between them should be 0.
  
  expect_equal(TquestDistance(x, y, tx, ty, tau=6), 0)
  
  # This value should not change if we change tx and ty. For example, 
  # if we do not provide any values, then they are assumed to be 
  #tx[i]={i/L, where L is the lenght of the series] 
  
  expect_equal(TquestDistance(x, y,  tau=6), 0)
  
  # Same if one of tx or tx is missing. In this case, since the length of the 
  # series is equal, the same sampling rate will be assigned.
  
  expect_equal(TquestDistance(x, y, tx=tx, tau=6), 0)
  expect_equal(TquestDistance(x, y, ty=ty, tau=6), 0)

  
  # Second example
  x <- c(5, 5, 7, 7, 7, 5, 5)
  y <- c(2, 8, 8, 8, 8, 2, 2)
  tx <- c(1, 2, 3, 4, 5, 6, 7)
  ty <- c(1, 2, 3, 4, 5, 6, 7)
  
  # Now, x trespasses tau in intervals {3,4,5} but y trespasses tau in 
  # {2,3,4,5}. Following the equation for the TQuest distance the distance between 
  # them is 1/1 * 1+ 1/1 * 1 = 2
  
  expect_equal(TquestDistance(x, y, tx, ty, tau=6), 2)
  
  # Third example: if both series never trespass the theshold, the distance is 0
  
  x <- c(5, 5, 7, 7, 7, 5, 5)
  y <- c(2, 8, 8, 8, 8, 2, 2)

  expect_equal(TquestDistance(x, y, tau=11), 0)

  # Fourth example: if one of the series never trespasses the 
  #threshold and the other does, then the distance between them is infinity.
  
  x <- c(5, 5, 12, 12, 12, 5, 5)
  y <- c(2, 8, 8, 8, 8, 2, 2)

  expect_equal(TquestDistance(x, y, tau=11), Inf)
  
  x <- c(5, 5, 7, 7, 7, 5, 5)
  y <- c(2, 12, 12, 12, 12, 2, 2)

  expect_equal(TquestDistance(x, y, tau=11), Inf)
  
  #Fifth example: if the last point of the series is included in one of 
  # the intervals that trespass the threshold. Following the equation of the TQuest functions we have  1/1 * 1+ 1/1 * 1 = 2.
  
  x <- c(5, 5, 7, 7, 7, 7, 7)
  y <- c(2, 12, 12, 12, 12, 12, 12)
  tx <- c(1, 2, 3, 4, 5, 6, 7)
  ty <- c(1, 2, 3, 4, 5, 6, 7)
  expect_equal(TquestDistance(x, y, tx, ty, tau=6), 2)
  
  y <- c(5, 5, 7, 7, 7, 7, 7)
  x <- c(2, 12, 12, 12, 12, 12, 12)
  ty <- c(1, 2, 3, 4, 5, 6, 7)
  tx <- c(1, 2, 3, 4, 5, 6, 7)
  expect_equal(TquestDistance(x, y, tx, ty, tau=6), 2)


  #Sixth example: case where x and y go over the threshold in more than
  #one interval. Following the equation of the TQuest functions we have  
  #1/2 * 1+ 1/2 * 1 = 1.
  
  x <- c(5, 5, 7, 7, 2, 2, 8, 8)
  y <- c(4, 4, 9, 9, 9, 1, 10, 10)
  tx <- c(1, 2, 3, 4, 5, 6, 7, 8)
  ty <- c(1, 2, 3, 4, 5, 6, 7, 8)
  expect_equal(TquestDistance(x, y, tx, ty, tau=6), 1)

})


test_that("Exceptions in TQuest distance", {
  x <- c("a","b","c","d")
  y <- c(3, 4, 1, 2)
  tau <- 2
  expect_equal(TquestDistance(x, y, tau=tau), NA)
  
  x <- replicate(3, rnorm(3)) 
  expect_equal(TquestDistance(x, y, tau=tau), NA)
  
  x <- c(1)
  expect_equal(TquestDistance(x, y, tau=tau), NA)
   
  x <- c(1, 2, NA, 3)
  expect_equal(TquestDistance(x, y, tau=tau), NA)
  
  x <- c(1, 2, 3, 4)
  tau <- "a"
  expect_equal(TquestDistance(x, y, tau=tau), NA)
  
  tx <- c(1, 2, 3, 4)
  ty <- c(2, 3, 4, NA)
  tau <- 1
  expect_equal(TquestDistance(x, y, tx, ty, tau=tau), NA)

  tx <- c(1, 2, 3, 4)
  ty <- c(2, 3, 4)
  expect_equal(TquestDistance(x, y, tx, ty, tau=tau), NA)
  
  tx <- c(1, -2, 3, 4)
  ty <- c(1, 2, 3, 4)
  expect_equal(TquestDistance(x, y, tx, ty, tau=tau), NA)
  
  tx <- c(1, 3, 2, 4)
  expect_equal(TquestDistance(x, y, tx, ty, tau=tau), NA)
  
  
})
