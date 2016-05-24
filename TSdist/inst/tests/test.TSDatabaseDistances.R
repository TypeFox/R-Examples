context("TSDatabaseDistance")

test_that("The TSDatabaseDistance function is checked", {
  
  x <- cumsum(rnorm(50))
  X <- rbind(x, x)
  
  #  Calculation of all distances for two identical series should be 0, except:
  
    # LCSS which is equal to the length of the series.
    # CMD: The distance betweem two equal series is not 0 (see reference)
    # NCD: Similar to the previous case, 
    # Pred: The distance is based on forecast densities, which depend on some random variable.
  
  
  # # # # # # # # # # # # # # # # # # # # # # # # FOR ONLY ONE DATABASE# # # # # # # # # # # # # # # # # # # # # # 
  
  # X must be a matrix, zoo, mts, xts or list
  expect_error(TSDatabaseDistances(X="a", distance="euclidean"))
  

  # If X is given as a ts object
  X1 <- as.ts(t(X))
  expect_equal(as.numeric(TSDatabaseDistances(X1, "euclidean")), 0)

  # It must contain more than one series
  X1 <- as.matrix(t(x))
  expect_error(TSDatabaseDistances(X1, "euclidean"))
  
  #  If X is given as a zoo object
  X1 <- as.zoo(t(X))
  expect_equal(as.numeric(TSDatabaseDistances(X1, "euclidean")), 0)
  
  # If X is given as a xts object
  data(zoo.database)
  X1 <- as.xts(cbind(zoo.database[,1], zoo.database[,1])) 
  expect_equal(as.numeric(TSDatabaseDistances(X1, "euclidean")), 0)
  
  # If X is given as a list
  X1 <- list(x, x)
  expect_equal(as.numeric(TSDatabaseDistances(X1, "euclidean")), 0)
  
  # If the length of the list is only one, an error is issued
  X1 <- list(x)
  expect_error(TSDatabaseDistances(X1, "euclidean"))

  # We check all distances with only one dataset
  expect_equal(as.numeric(TSDatabaseDistances(X, "euclidean")), 0)
  expect_equal(as.numeric(TSDatabaseDistances(X, "manhattan")), 0)
  expect_equal(as.numeric(TSDatabaseDistances(X, "minkowski", p=3)), 0)
  expect_equal(as.numeric(TSDatabaseDistances(X, "infnorm")), 0)
  expect_equal(as.numeric(TSDatabaseDistances(X, "ccor")), 0)
  expect_equal(as.numeric(TSDatabaseDistances(X, "sts")), 0)
  expect_equal(as.numeric(TSDatabaseDistances(X, "dtw")), 0)
  expect_equal(as.numeric(TSDatabaseDistances(X, "lb.keogh", window.size=3)), 0)
  expect_equal(as.numeric(TSDatabaseDistances(X, "edr", epsilon=0.1)), 0)
  expect_equal(as.numeric(TSDatabaseDistances(X, "erp", g=0.1)), 0)
  expect_equal(as.numeric(TSDatabaseDistances(X, "lcss", epsilon=0.1)), 50)
  expect_equal(as.numeric(TSDatabaseDistances(X, "fourier")), 0)
  expect_equal(as.numeric(TSDatabaseDistances(X, "tquest", tau=0.1)), 0)
  expect_equal(as.numeric(TSDatabaseDistances(X, "dissim")), 0)
  expect_equal(as.numeric(TSDatabaseDistances(X, "acf")), 0)
  expect_equal(as.numeric(TSDatabaseDistances(X, "pacf")), 0)
  expect_equal(as.numeric(TSDatabaseDistances(X, "ar.lpc.ceps")), 0)
  expect_equal(as.numeric(as.numeric(TSDatabaseDistances(X, "ar.mah")$statistic)), 0)
  expect_equal(as.numeric(TSDatabaseDistances(X, "ar.pic")), 0) 
  expect_equal(as.numeric(TSDatabaseDistances(X, "cdm")), CDMDistance(x,x))
  expect_equal(as.numeric(TSDatabaseDistances(X, "cid")), 0)
  expect_equal(round(as.numeric(TSDatabaseDistances(X, "cor")), digits=2), 0)
  expect_equal(as.numeric(TSDatabaseDistances(X, "cort")), 0)
  expect_equal(as.numeric(TSDatabaseDistances(X, "wav")), 0)
  expect_equal(as.numeric(TSDatabaseDistances(X, "int.per")), 0)
  expect_equal(as.numeric(TSDatabaseDistances(X, "per")), 0)
  expect_equal(as.numeric(TSDatabaseDistances(X, "mindist.sax", w=2)), 0)
  expect_equal(as.numeric(TSDatabaseDistances(X, "ncd")), NCDDistance(x,x))
  set.seed(123)
  d <- PredDistance(x,x, h=3)
  set.seed(123)
  expect_equal(as.numeric(TSDatabaseDistances(X, "pred", h=3)), d)
  expect_equal(as.numeric(TSDatabaseDistances(X, "spec.glk")), 0)
  expect_equal(as.numeric(TSDatabaseDistances(X, "spec.isd")), 0)
  expect_equal(as.numeric(TSDatabaseDistances(X, "spec.llr")), 0)
  expect_equal(as.numeric(TSDatabaseDistances(X, "pdc")), 0)
  expect_equal(as.numeric(TSDatabaseDistances(X, "frechet")), 0)
  
  # Some other options that must be explored

  expect_equal(as.numeric(TSDatabaseDistances(X, "spec.llr", method="LK")), 0)
  expect_equal(as.numeric(TSDatabaseDistances(X, "spec.llr", n=0)), 0)
  expect_equal(as.numeric(TSDatabaseDistances(X, "spec.isd", n=0)), 0)
  
  # If X is a list a different treatment is necessary for some distances
  X1 <- list(x, x)
  expect_equal(as.numeric(TSDatabaseDistances(X1, "ar.pic")), 0)
  
  # Other exceptions specific to some distance measures
  order.x <- matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0), nrow=3)
  expect_error(as.numeric(TSDatabaseDistances(X, "ar.pic", order.x=order.x)))
  expect_error(as.numeric(TSDatabaseDistances(X, "ar.lpc.ceps", order.x=order.x)))
  
  logarithms.x <- matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0), nrow=3)
  expect_error(as.numeric(TSDatabaseDistances(X, "pred", logarithms.x=logarithms.x)))
  
  differences.x <- matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0), nrow=3)
  expect_error(as.numeric(TSDatabaseDistances(X, Y, "pred", differences.x=differences.x)))
  
  expect_error(as.numeric(TSDatabaseDistances(X, "spec.llr", method="AA")))
  
  #  An special if, for when the number of series is uneven 
  X <- rbind(x, x, x)
  set.seed(123)
  d <- as.numeric(diss(X, METHOD="PRED",h=3)$dist)
  set.seed(123)
  d1 <-as.numeric(t(TSDatabaseDistances(X, "pred", h=3)))
  expect_equal(d1, d)
  
  # # # # # # # # # # # # # # # # # # FOR TWO DATASETS# # # # # # # # # # # # # # # # # # # # # # # # # # # 
  
  X <- rbind(x, x)
  Y <- X
  y <- cumsum(rnorm(50))
  M <- c(0, 0, 0, 0)
  
  # X and Y must be a matrix, zoo, mts, xts or list
  expect_error(TSDatabaseDistances(X=X, Y=mean, distance="euclidean"))
  
  # If Y is given as a ts object
  Y1 <- as.ts(t(Y))
  expect_equal(as.numeric(TSDatabaseDistances(X, Y1, "euclidean")), M)
  
  # It must contain more than one series
  Y1 <- as.matrix(t(y))
  expect_error(TSDatabaseDistances(X, Y1, "euclidean"))
  
  # If Y is given as a zoo object
  Y1 <- as.zoo(t(Y))
  expect_equal(as.numeric(TSDatabaseDistances(X, Y1, "euclidean")), M)
  
  # If Y is given as a xts object
  data(zoo.database)
  X1 <- as.xts(cbind(zoo.database[,1], zoo.database[,1])) 
  Y1 <- as.xts(cbind(zoo.database[,1], zoo.database[,1])) 
  expect_equal(as.numeric(TSDatabaseDistances(X1, Y1, "euclidean")), M)
  
  # If Y is given as a list
  Y1 <- list(x, x)
  expect_equal(as.numeric(TSDatabaseDistances(Y1, "euclidean")), 0)
  
  # If the length of the list is only one, an error is issued
  Y1 <- list(x)
  expect_error(TSDatabaseDistances(X, Y1, "euclidean"))
  
  # We calculate all the distance measures for two matrices  
  expect_equal(as.numeric(TSDatabaseDistances(X, Y, "euclidean")), M)
  expect_equal(as.numeric(TSDatabaseDistances(X, Y, "manhattan")), M)
  expect_equal(as.numeric(TSDatabaseDistances(X, Y, "minkowski", p=3)), M)
  expect_equal(as.numeric(TSDatabaseDistances(X, Y, "infnorm")), M)
  expect_equal(as.numeric(TSDatabaseDistances(X, Y, "ccor")), M)
  expect_equal(as.numeric(TSDatabaseDistances(X, Y, "sts")), M)
  expect_equal(as.numeric(TSDatabaseDistances(X, Y, "dtw")), M)
  expect_equal(as.numeric(TSDatabaseDistances(X, Y, "lb.keogh", window.size=3)), M)
  expect_equal(as.numeric(TSDatabaseDistances(X, Y, "edr", epsilon=0.1)), M)
  expect_equal(as.numeric(TSDatabaseDistances(X, Y, "erp", g=0.1)), M)
  expect_equal(as.numeric(TSDatabaseDistances(X, Y, "lcss", epsilon=0.1)), c(50, 50, 50, 50))
  expect_equal(as.numeric(TSDatabaseDistances(X, Y, "fourier")), M)
  expect_equal(as.numeric(TSDatabaseDistances(X, Y, "tquest", tau=0.1)), M)
  expect_equal(as.numeric(TSDatabaseDistances(X, Y, "dissim")), M)
  expect_equal(as.numeric(TSDatabaseDistances(X, Y, "acf")), M)
  expect_equal(as.numeric(TSDatabaseDistances(X, Y, "pacf")), M)
  expect_equal(as.numeric(TSDatabaseDistances(X, Y, "ar.lpc.ceps")), M)
  expect_equal(as.numeric(as.numeric(TSDatabaseDistances(X, Y, "ar.mah")$statistic)), M)
  expect_equal(as.numeric(TSDatabaseDistances(X, Y, "ar.pic")), M)
  expect_equal(as.numeric(TSDatabaseDistances(X, Y, "cdm")), rep(CDMDistance(x,x),4))  
  expect_equal(as.numeric(TSDatabaseDistances(X, Y, "cid")), M)
  expect_equal(round(as.numeric(TSDatabaseDistances(X, Y, "cor")), digits=2), M)
  expect_equal(as.numeric(TSDatabaseDistances(X, Y, "cort")), M)
  expect_equal(as.numeric(TSDatabaseDistances(X, Y, "wav")), M)
  expect_equal(as.numeric(TSDatabaseDistances(X, Y, "int.per")), M)
  expect_equal(as.numeric(TSDatabaseDistances(X, Y, "per")), M)
  expect_equal(as.numeric(TSDatabaseDistances(X, Y, "mindist.sax", w=2)), M)
  expect_equal(as.numeric(TSDatabaseDistances(X, Y, "ncd")), rep(NCDDistance(x,x), 4))
  set.seed(123)
  d <- diss(SERIES=rbind(X,X), METHOD="PRED", h=3)$dist
  d <- as.numeric(d)
  d <- d[c(2, 4, 3, 5)]
  set.seed(123)
  expect_equal(as.numeric(TSDatabaseDistances(X, X, "pred", h=3)), d)
  expect_equal(as.numeric(TSDatabaseDistances(X, Y, "spec.glk")), M)
  expect_equal(as.numeric(TSDatabaseDistances(X, Y, "spec.isd")), M)
  expect_equal(as.numeric(TSDatabaseDistances(X, Y, "spec.llr")), M)
  expect_equal(as.numeric(TSDatabaseDistances(X, Y, "pdc")), M)
  expect_equal(as.numeric(TSDatabaseDistances(X, X, "frechet")), M)
  
  # If X is a list a different treatment is necessary for some distances
  X1 <- list(x, x)
  Y1 <- list(x, x)
  expect_equal(as.numeric(TSDatabaseDistances(X1, Y1, "wav")), M)
  expect_equal(as.numeric(TSDatabaseDistances(X1, Y1, "ar.pic")), M)
  
  # Other exceptions specific to some distance measures
  order.y <- matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0), nrow=3)
  expect_error(as.numeric(TSDatabaseDistances(X, Y, "ar.pic", order.y=order.y)))
  expect_error(as.numeric(TSDatabaseDistances(X, Y, "ar.lpc.ceps", order.y=order.y)))
  
  logarithms.x <- matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0), nrow=3)
  logarithms.y <- matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0), nrow=3)
  expect_error(as.numeric(TSDatabaseDistances(X, Y, "pred", logarithms.x=logarithms.x)))
  expect_error(as.numeric(TSDatabaseDistances(X, Y, "pred", h=3, logarithms.y=logarithms.y)))
  
  differences.x <- matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0), nrow=3)
  differences.y <- matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0), nrow=3)
  expect_error(as.numeric(TSDatabaseDistances(X, Y, "pred", differences.x=differences.x)))
  expect_error(as.numeric(TSDatabaseDistances(X, Y, "pred", h=3, differences.y=differences.y)))
  
  expect_error(as.numeric(TSDatabaseDistances(X, X, "pred")))
  
  #  Some other special cases must be explored
  expect_equal(as.numeric(TSDatabaseDistances(X, Y, "spec.llr", n=0)), rep(0, 4))
  
  expect_equal(as.numeric(TSDatabaseDistances(X, Y, "spec.isd", n=0)), rep(0, 4))
  
  #  An special if, for when the number of series is uneven 
  X <- rbind(x, x, x)
  Y <- rbind(x, x)
  set.seed(123)
  d <- diss(rbind(X,Y), METHOD="PRED",h=3)$dist
  d <- d[c(3, 4, 6, 7, 8, 9)]
  set.seed(123)
  d1 <-as.numeric(t(TSDatabaseDistances(X, Y, "pred", h=3)))
  expect_equal(d1, d)
  

  
  
})
