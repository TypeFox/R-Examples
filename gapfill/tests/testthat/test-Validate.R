#require(testthat); library("gapfill", lib.loc = "../../../lib/")
context("test-Validate")

test_that("Validate",{

  dataO <- c(1:10, rep(NA, 10))
  dataP <- c(rep(NA, 10), 2:11)
  dataT <- c(1:10, 1:10)
  t1 <- cbind(nNA = 10,
              nFilled = 10,
              nNotFilled = 0,
              ratioFilled = 1,
              nCrossvali = 10,
              RMSE = 1)
  
  expect_equal(Validate(dataO, dataP, dataT), t1)
  expect_equal(Validate(dataO, dataP, dataT,
                        include = rep(c(FALSE, TRUE), c(10, 10))),
               t1)
  t2 <- cbind(nNA = 10,
              nFilled = 10,
              nNotFilled = 0,
              ratioFilled = 1,
              nCrossvali = 5,
              RMSE = 1)
  expect_equal(Validate(dataO, dataP, dataT,
                        include = rep(c(FALSE, TRUE), c(15, 5))),
               t2)

  dataP <- c(rep(NA, 12), 4:11)
  t3 <- cbind(nNA = 10,
              nFilled = 8,
              nNotFilled = 2,
              ratioFilled = 8/10,
              nCrossvali = 8,
              RMSE = 1)
  expect_equal(Validate(dataO, dataP, dataT), t3)

  t3 <- cbind(nNA = 10,
              nFilled = 8,
              nNotFilled = 2,
              ratioFilled = 8/10,
              nCrossvali = 5,
              RMSE = 1)
  expect_equal(Validate(dataO, dataP, dataT,
                        include = rep(c(FALSE, TRUE), c(15, 5))),
               t3)
  
  vp <- 300 + c(5:10) + rep(21 * c(0:5), each = 6)
  nn <- ndvi
  nn[vp] <- NA
  out <- Gapfill(nn, subset = vp)
  expect_equal(Validate(dataObserved = nn,  dataFilled = out$fill,
                        dataTrue = ndvi),
               structure(c(1639, 36, 1603, 0.0219646125686394, 36, 0.0222318373858305
), .Dim = c(1L, 6L), .Dimnames = list(NULL, c("nNA", "nFilled", 
"nNotFilled", "ratioFilled", "nCrossvali", "RMSE"))),
               tolerance = 1e-4)
})

