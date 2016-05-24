context("TSClust")

test_that("The ACF and PACF distance wrappers work properly", {
  
  # We use examples in the documentation of the TSclust package to test 
  # that the wrapper works correctly: 
  
  x <- cumsum(rnorm(100))
  y <- cumsum(rnorm(100))
  z <- sin(seq(0, pi, length.out=100))
  
  d1 <- as.numeric(diss.PACF(x, y))
  d2 <- as.numeric(diss.ACF(x, z))
  d3 <- as.numeric(diss.PACF(y, z))
  
  expect_equal(PACFDistance(x, y), d1)
  expect_equal(ACFDistance(x, z), d2)
  expect_equal(PACFDistance(y, z), d3)
  
  # If an error is thrown by the original function, the wrapper will return NA.
  # For example, if there are missing values in the series:
  x[1] <- NA
  expect_equal(ACFDistance(x, y), NA)
  expect_equal(PACFDistance(x, y), NA)
  
  # For specific errors of the original function access the documentation of the
  # TSclust package.

})


test_that("The AR LPC CEPS distance wrapper works properly", {
  
  # We use examples in the documentation of the TSclust package to test that the
  # wrapper works correctly: 

  x <- arima.sim(model=list(ar=c(0.4,-0.1)), n =100, n.start=100)
  y <- arima.sim(model=list(ar=c(0.9)), n =100, n.start=100)
  z <- arima.sim(model=list(ar=c(0.5, 0.2)), n =100, n.start=100)
  
  d1 <- diss.AR.LPC.CEPS(x, y, 25)
  d2 <- diss.AR.LPC.CEPS(x, z, 25, order.x = c(2,0,0), order.y = NULL )
  d3 <- diss.AR.LPC.CEPS(y, z, 25)
  
  expect_equal(ARLPCCepsDistance(x, y, 25), d1)
  expect_equal(ARLPCCepsDistance(x, z, 25, order.x = c(2,0,0), order.y = NULL ), 
               d2)
  expect_equal(ARLPCCepsDistance(y, z, 25), d3)
  
  # If an error is thrown by the original function, the wrapper will return NA.
  # For example, if there are missing values in the series:
  x[1] <- NA
  expect_equal(ARLPCCepsDistance(x, y, 25), NA)
  
  # For specific errors of the original function access the documentation of the
  # TSclust package.

})


test_that("The AR MAH distance wrapper works properly", {
  
  # We use examples in the documentation of the TSclust package to test that the
  # wrapper works correctly: 

  x <- arima.sim(model=list(ar=c(0.4,-0.1)), n =100, n.start=100)
  y <- arima.sim(model=list(ar=c(0.9)), n =100, n.start=100)
  z <- arima.sim(model=list(ar=c(0.5, 0.2)), n =100, n.start=100)

  d1 <- diss.AR.MAH(x, y)
  d2 <- diss.AR.MAH(x, z)
  d3 <- diss.AR.MAH(y, z)
  
  expect_equal(ARMahDistance(x, y), d1)
  expect_equal(ARMahDistance(x, z), d2)
  expect_equal(ARMahDistance(y, z), d3)
  
  #For specific statistic and p-value distances
  expect_equal(ARMahStatisticDistance(x, y), digits=2, as.numeric(d1$statistic, digits=2))
  expect_equal(ARMahPvalueDistance(x, y), as.numeric(d1$p_value))
  
  # If an error is thrown by the original function, the wrapper will return NA.
  # For example, if there are missing values in the series:
  x[1] <- NA
  expect_equal(ARMahDistance(x, y), NA)
  expect_equal(ARMahStatisticDistance(x, y), NA)
  expect_equal(ARMahPvalueDistance(x, y), NA)
  
  # For specific errors of the original function access the documentation of the
  # TSclust package.

})


test_that("The AR PIC distance wrapper works properly", {
  
  # We use examples in the documentation of the TSclust package to test that the
  # wrapper works correctly: 

  x <- arima.sim(model=list(ar=c(0.4,-0.1)), n =100, n.start=100)
  y <- arima.sim(model=list(ar=c(0.9)), n =100, n.start=100)
  z <- arima.sim(model=list(ar=c(0.5, 0.2)), n =100, n.start=100)

  d1 <- diss.AR.PIC( x, y, order.x = c(2,0,0), order.y = c(1,0,0) )
  d2 <- diss.AR.PIC( x, z, order.x = c(2,0,0), order.y = c(2,0,0) )
  d3 <- diss.AR.PIC( y, z, order.x=NULL, order.y=c(2,0,0) )
  
  expect_equal(ARPicDistance(x, y, order.x = c(2,0,0), order.y = c(1,0,0)), d1)
  expect_equal(ARPicDistance(x, z, order.x = c(2,0,0), order.y = c(2,0,0)), d2)
  expect_equal(ARPicDistance(y, z, order.x=NULL, order.y=c(2,0,0) ), d3)
  
  # If an error is thrown by the original function, the wrapper will return NA.
  # For example, if there are missing values in the series:
  x[1] <- NA
  expect_equal(ARPicDistance(x, y, order.x = c(2,0,0), order.y = c(1,0,0)), NA)
  
  # For specific errors of the original function access the documentation of the
  # TSclust package.

})

test_that("The CDM distance wrapper works properly", {
  
  # We use examples in the documentation of the TSclust package to test that the
  # wrapper works correctly: 

  n = 50
  x <- rnorm(n) 
  y <- cumsum(rnorm(n))

  d1 <- diss.CDM(x, y)
  d2 <- diss.CDM(x, y, type="bzip2")
  
  expect_equal(CDMDistance(x, y), d1)
  expect_equal(CDMDistance(x, y, type="bzip2"), d2)
  
  # If an error is thrown by the original function, the wrapper will return NA.
  # For example, if there are missing values in the series:
  x[1] <- NA
  expect_equal(CDMDistance(x, y, type="bzip2"), NA)
  
  # For specific errors of the original function access the documentation of the
  # TSclust package.

})

test_that("The CID distance wrapper works properly", {
  
  # We use examples in the documentation of the TSclust package to test that the
  # wrapper works correctly: 

  n <- 50
  x <- rnorm(n) 
  y <- cumsum(rnorm(n))

  d1 <- as.numeric(diss.CID(x, y))
  
  expect_equal(CIDDistance(x, y), d1)
  
  # If an error is thrown by the original function, the wrapper will return NA.
  # For example, if there are missing values in the series:
  x[1] <- NA
  expect_equal(CIDDistance(x, y), NA)
  
  # For specific errors of the original function access the documentation of the
  # TSclust package.

})

test_that("The Cor distance wrapper works properly", {
  
  # We use examples in the documentation of the TSclust package to test that 
  # the wrapper works correctly: 
  
  x <- cumsum(rnorm(100))
  y <- cumsum(rnorm(100))
  z <- sin(seq(0, pi, length.out=100))
  
  d1 <- diss.COR(x, y)
  d2 <- diss.COR(x, z)
  d3 <- diss.COR(y, z)
  
  expect_equal(CorDistance(x, y), d1)
  expect_equal(CorDistance(x, z), d2)
  expect_equal(CorDistance(y, z), d3)
  
  # If an error is thrown by the original function, the wrapper will return NA.
  # For example, if there are missing values in the series:
  x[1] <- NA
  expect_equal(CorDistance(x, y), NA)
  
  # For specific errors of the original function access the documentation of the
  # TSclust package.
  
})


test_that("The Cort distance wrapper works properly", {
  
  # We use examples in the documentation of the TSclust package to test that the
  # wrapper works correctly: 

  x <- cumsum(rnorm(100))
  y <- cumsum(rnorm(100))
  z <- sin(seq(0, pi, length.out=100))

  d1 <- diss.CORT(x, y, 2)
  d2 <- diss.CORT(x, z, 2)
  d3 <- diss.CORT(y, z, 2)

  expect_equal(CortDistance(x, y), d1)
  expect_equal(CortDistance(x, z), d2)
  expect_equal(CortDistance(y, z), d3)
  
  # If an error is thrown by the original function, the wrapper will return NA.
  # For example, if there are missing values in the series:
  x[1] <- NA
  expect_equal(CortDistance(x, y), NA)
  
  # For specific errors of the original function access the documentation of the
  # TSclust package.

})

test_that("The Wavelet distance wrapper works properly", {
  
  # We use examples in the documentation of the TSclust package to test that the
  # wrapper works correctly: 

  x <- cumsum(rnorm(100))
  y <- cumsum(rnorm(100))
  z <- sin(seq(0, pi, length.out=100))

  d1 <- as.numeric(diss.DWT(rbind(x, y)))
  d2 <- as.numeric(diss.DWT(rbind(x, z)))
  d3 <- as.numeric(diss.DWT(rbind(y, z)))

  expect_equal(WavDistance(x, y), d1)
  expect_equal(WavDistance(x, z), d2)
  expect_equal(WavDistance(y, z), d3)
  
  # If an error is thrown by the original function, the wrapper will return NA.
  # For example, if there are missing values in the series:
  x[1] <- NA
  expect_equal(WavDistance(x, y), NA)
  
  # For specific errors of the original function access the documentation of the
  # TSclust package.

})


test_that("The Int Per distance wrapper works properly", {
  
  # We use examples in the documentation of the TSclust package to test that the
  # wrapper works correctly: 

  x <- cumsum(rnorm(100))
  y <- cumsum(rnorm(100))
  z <- sin(seq(0, pi, length.out=100))

  d1 <- diss.INT.PER(x, y, normalize=TRUE)
  d2 <- diss.INT.PER(x, z, normalize=TRUE)
  d3 <- diss.INT.PER(y, z, normalize=TRUE)
  
  expect_equal(IntPerDistance(x, y, normalize=TRUE), d1)
  expect_equal(IntPerDistance(x, z, normalize=TRUE), d2)
  expect_equal(IntPerDistance(y, z, normalize=TRUE), d3)
  
  # If an error is thrown by the original function, the wrapper will return NA.
  # For example, if there are missing values in the series:
  x[1] <- NA
  expect_equal(IntPerDistance(x, y, normalize=TRUE), NA)
  
  # For specific errors of the original function access the documentation of the
  # TSclust package.

})


test_that("The Per distance wrapper works properly", {
  
  # We use examples in the documentation of the TSclust package to test that the
  # wrapper works correctly: 

  x <- cumsum(rnorm(100))
  y <- cumsum(rnorm(100))
  z <- sin(seq(0, pi, length.out=100))

  d1 <- as.numeric(diss.PER(x, y))
  d2 <- as.numeric(diss.PER(x, z))
  d3 <- as.numeric(diss.PER(y, z))
  d4 <- as.numeric(diss.PER(x, y, TRUE, TRUE))
  d5 <- as.numeric(diss.PER(x, z, TRUE, TRUE))
  d6 <- as.numeric(diss.PER(y, z, TRUE, TRUE))
  
  expect_equal(PerDistance(x, y), d1)
  expect_equal(PerDistance(x, z), d2)
  expect_equal(PerDistance(y, z), d3)
  expect_equal(PerDistance(x, y, TRUE, TRUE), d4)
  expect_equal(PerDistance(x, z, TRUE, TRUE), d5)
  expect_equal(PerDistance(y, z, TRUE, TRUE), d6)
  
  # If an error is thrown by the original function, the wrapper will return NA.
  # For example, if there are missing values in the series:
  x[1] <- NA
  expect_equal(PerDistance(x, y), NA)
  
  # For specific errors of the original function access the documentation of the
  # TSclust package.

})


test_that("The Mindist SAX distance wrapper works properly", {
  
  # We use examples in the documentation of the TSclust package to test that the
  # wrapper works correctly: 

  set.seed(12349)
  n <- 100
  x <- rnorm(n)  #generate sample series, white noise and a wiener process
  y <- cumsum(rnorm(n))
  z <- rnorm(n)^2
  w <- 20 #amount of equal-sized frames to divide the series
  alpha <- 4 #size of the alphabet, parameter for 
  

  d1 <- diss.MINDIST.SAX(x, y, w, alpha, plot=FALSE)
  d2 <- diss.MINDIST.SAX(x, z, w, alpha, plot=FALSE)
  d3 <- diss.MINDIST.SAX(y, z, w, alpha, plot=FALSE)

  
  expect_equal(MindistSaxDistance(x, y, w, alpha, plot=FALSE), d1)
  expect_equal(MindistSaxDistance(x, z, w, alpha, plot=FALSE), d2)
  expect_equal(MindistSaxDistance(y, z, w, alpha, plot=FALSE), d3)
  
  # If an error is thrown by the original function, the wrapper will return NA.
  # For example, if there are missing values in the series:
  x[1] <- NA
  expect_equal(MindistSaxDistance(x, y, w, alpha, plot=FALSE), NA)
  
  # For specific errors of the original function access the documentation of the
  # TSclust package.

})


test_that("The NCD distance wrapper works properly", {
  
  # We use examples in the documentation of the TSclust package to test that the
  # wrapper works correctly: 

  n = 50
  x <- rnorm(n)  
  y <- cumsum(rnorm(n))

  d1 <- diss.NCD(x, y)
  d2 <- diss.NCD(x, y, type="bzip2")
  
  expect_equal(NCDDistance(x, y), d1)
  expect_equal(NCDDistance(x, y, type="bzip2"), d2)
  
  # If an error is thrown by the original function, the wrapper will return NA.
  # For example, if there are missing values in the series:
  x[1] <- NA
  expect_equal(NCDDistance(x, y, type="bzip2"), NA)
  
  # For specific errors of the original function access the documentation of the
  # TSclust package.

})


test_that("The Pred distance wrapper works properly", {
  
  # We use examples in the documentation of the TSclust package to test that the
  # wrapper works correctly: 

  x <- (rnorm(100))
  x <- x + abs(min(x)) + 1 
  y <- (rnorm(100))

  set.seed(123)
  d1 <- diss.PRED(x, y, h=6, logarithm.x=FALSE, logarithm.y=FALSE, 
                  differences.x=1, differences.y=0)$L1dist
  
  set.seed(123)
  expect_equal(PredDistance(x, y, h=6, logarithm.x=FALSE, logarithm.y=FALSE, 
                            differences.x=1, differences.y=0), d1)
  
  # If an error is thrown by the original function, the wrapper will return NA.
  # For example, if the length of differences.x is not equal to
  # the number of the series in the database.
  x[1] <- NA
  expect_equal(PredDistance(x, y,  h=6, logarithm.x=FALSE, logarithm.y=FALSE, differences.x=1, differences.y=0), NA)
  
  # For specific errors of the original function access the documentation of the
  # TSclust package.

})

test_that("The Spec GLK distance wrapper works properly", {
  
  # We use examples in the documentation of the TSclust package to test that the
  # wrapper works correctly: 

  x <- cumsum(rnorm(50))
  y <- cumsum(rnorm(50))
  z <- sin(seq(0, pi, length.out=50))
  
  #In this case we use diss because diss.SPEC.GLK has a bug and always returns 0.
  d1 <- as.numeric(diss(rbind(x, y), "SPEC.GLK"))
  d2 <- as.numeric(diss(rbind(x, z), "SPEC.GLK"))
  d3 <- as.numeric(diss(rbind(y, z), "SPEC.GLK"))

  expect_equal(SpecGLKDistance(x, y), d1)
  expect_equal(SpecGLKDistance(x, z), d2)
  expect_equal(SpecGLKDistance(y, z), d3)
  
  # If an error is thrown by the original function, the wrapper will return NA.
  # For example, if there are missing values in the series:
  x[1] <- NA
  expect_equal(SpecGLKDistance(x, y), NA)
  
  # For specific errors of the original function access the documentation of the
  # TSclust package.

})

test_that("The Spec ISD distance wrapper works properly", {
  
  # We use examples in the documentation of the TSclust package to test that the
  # wrapper works correctly: 

  x <- cumsum(rnorm(50))
  y <- cumsum(rnorm(50))
  z <- sin(seq(0, pi, length.out=50))
  
  d1 <- as.numeric(diss.SPEC.ISD(x, y))
  d2 <- as.numeric(diss.SPEC.ISD(x, z))
  d3 <- as.numeric(diss.SPEC.ISD(y, z))

  expect_equal(SpecISDDistance(x, y), d1)
  expect_equal(SpecISDDistance(x, z), d2)
  expect_equal(SpecISDDistance(y, z), d3)

  # If an error is thrown by the original function, the wrapper will return NA.
  # For example, if there are missing values in the series:
  x[1] <- NA
  expect_equal(SpecISDDistance(x, y), NA)
  
  # For specific errors of the original function access the documentation of the
  # TSclust package.
  
})

test_that("The Spec LLR distance wrapper works properly", {
  
  # We use examples in the documentation of the TSclust package to test that the
  # wrapper works correctly: 

  x <- cumsum(rnorm(50))
  y <- cumsum(rnorm(50))
  z <- sin(seq(0, pi, length.out=50))
  

  d1 <- as.numeric(diss.SPEC.LLR(x, y))
  d2 <- as.numeric(diss.SPEC.LLR(x, z))
  d3 <- as.numeric(diss.SPEC.LLR(y, z))

  expect_equal(SpecLLRDistance(x, y), d1)
  expect_equal(SpecLLRDistance(x, z), d2)
  expect_equal(SpecLLRDistance(y, z), d3)
  
  # If an error is thrown by the original function, the wrapper will return NA.
  # For example, if there are missing values in the series:
  x[1] <- NA
  expect_equal(SpecLLRDistance(x, y), NA)
  
  # For specific errors of the original function access the documentation of the
  # TSclust package.
  

})

