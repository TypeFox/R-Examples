library(dsm)
library(Distance)
library(testthat)

lnl.tol<-1e-4
par.tol<-1e-6

context("Mexico pantropical dolphin data")

# load the Gulf of Mexico dolphin data
data(mexdolphins)

# fit a detection function and look at the summary
hn.model <- suppressMessages(ds(mexdolphins$distdata,
                                max(mexdolphins$distdata$distance),
                                adjustment = NULL))

test_that("Do we get the same results?",{


  ddf.par <- 8.626542
  names(hn.model$ddf$par) <- NULL
  expect_that(hn.model$ddf$par, equals(ddf.par,tolerance=par.tol))

  # fit a simple smooth of x and y
  mod1<-dsm(N~s(x,y), hn.model, mexdolphins$segdata, mexdolphins$obsdata)
  #summary(mod1)


  expect_that(mod1$gcv.ubre, equals(42.9169051,tolerance=par.tol))


  expect_that(dsm.cor(mod1,resid.type="d",max.lag=9),
              throws_error("No column called Segment.Label in data"))
})

test_that("Density weighting",{

  # fit density model
  mod1 <- dsm(D~s(x,y), hn.model, mexdolphins$segdata, mexdolphins$obsdata)

  # compare when we set the weights
  mod1.w <- dsm(D~s(x,y), hn.model, mexdolphins$segdata, mexdolphins$obsdata,
                weights=mod1$data$segment.area)

  expect_equal(fitted(mod1),fitted(mod1.w),tolerance=par.tol)


  # setting weights to 1 or another constant
  # compare when we set the weights
  mod1.w1 <- dsm(D~s(x,y), hn.model, mexdolphins$segdata, mexdolphins$obsdata,
                weights=rep(1,nrow(mexdolphins$segdata)))
  # compare when we set the weights
  mod1.w2 <- dsm(D~s(x,y), hn.model, mexdolphins$segdata, mexdolphins$obsdata,
                weights=rep(100,nrow(mexdolphins$segdata)))

  expect_equal(fitted(mod1.w1),fitted(mod1.w2),tolerance=par.tol)

  # scalar input of weights (same as weighting all as 1, or 10)
  mod1.ws1 <- dsm(D~s(x,y), hn.model, mexdolphins$segdata, mexdolphins$obsdata,
                weights=1)

  expect_equal(fitted(mod1.ws1),fitted(mod1.w2),tolerance=par.tol)

})

