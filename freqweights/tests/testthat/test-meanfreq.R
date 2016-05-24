## Emilio Torres Manzanera
## University of Oviedo
## Time-stamp: <2014-04-24 15:26 emilio on emilio-Satellite-P100>
## ============================================================

context("Mean")

test_that("Mean",{
  tfq <- tablefreq(iris$Sepal.Length)
  
  m <- mean(iris$Sepal.Length)
  a <- meanfreq(iris$Sepal.Length)
  b <- meanfreq(tfq,freq="freq")
  c <- .meanfreq(tfq)

  expect_that(a, equals(m))
  expect_that(b, equals(m))
  expect_that(c, equals(m))

  if(require(hflights)){
    tfq <- tablefreq(hflights$ArrDelay)

   m <- mean(hflights$ArrDelay, na.rm=TRUE)
  a <- meanfreq(hflights$ArrDelay)
  b <- meanfreq(tfq,freq="freq")
  c <- .meanfreq(tfq)

  expect_that(a, equals(m))
  expect_that(b, equals(m))
  expect_that(c, equals(m))
 }

  
})
