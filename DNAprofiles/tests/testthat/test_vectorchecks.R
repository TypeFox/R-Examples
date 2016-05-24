library(DNAprofiles)
context("Vector checks")

test_that(desc = "Vector checks",{
  
  x7 <- x6 <- x5 <- x4 <- x3 <- x2 <- x1 <- x <- seq(1,10,length.out = 40) 
  
  expect_true(Zallfinitepos(x,1,length(x)))  
  
  x1[3] <- Inf    
  expect_false(Zallfinitepos(x1,1,length(x))) 

  x2[3] <- -Inf
  expect_false(Zallfinitepos(x2,1,length(x))  )
  expect_true(Zallfinitepos(x2,1,2)  )
  expect_true(Zallfinitepos(x2,4,length(x))  )
  
  x3[4] <- NaN
  expect_false(Zallfinitepos(x3,1,length(x))  )
  
  x4[length(x4)] <- NA  
  expect_false(Zallfinitepos(x4,1,length(x)))
  
  x5[1] <- NA  
  expect_false(Zallfinitepos(x5,1,length(x)))
  
  x6[5] <- 0
  expect_false(Zallfinitepos(x6,1,length(x)))  
  
  x7[1] <- .Machine$double.xmax
  expect_true(Zallfinitepos(x7,1,length(x)))    
})