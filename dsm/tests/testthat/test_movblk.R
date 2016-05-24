library(dsm)
library(Distance)
library(testthat)

cv.tol<-1e-5
N.tol<-1e-4


## NB the 444km^2 for the prediction grid is INCORRECT but
## serves us fine for the purpose of these tests.

context("Moving block bootstrap")

# load the Gulf of Mexico dolphin data
data(mexdolphins)

# fit a detection function and look at the summary
hn.model <- suppressMessages(ds(mexdolphins$distdata,
                                max(mexdolphins$distdata$distance),
                                adjustment = NULL))

# fit a simple smooth of x and y
mod1<-dsm(N~s(x,y), hn.model, mexdolphins$segdata, mexdolphins$obsdata)

# run the moving block bootstrap for 2 rounds
set.seed(1123)
mod1.movblk <- dsm.var.movblk(mod1, mexdolphins$preddata, n.boot = 2,
                              block.size = 3, off.set = 444*1000*1000,bar=FALSE)

test_that("mexdolphins - bootstrap results for s(x,y)",{

  expect_that(mod1.movblk$study.area.total,
              equals(c(37065.62, 20722.39),tol=N.tol))

  expect_that(summary(mod1.movblk)$cv[1],
              equals(0.3367496,tol=cv.tol))

  expect_that(summary(mod1.movblk)$bootstrap.cv[1],
              equals(0.3117824,tol=cv.tol))
})

## With no detection function
test_that("mexdolphins - bootstrap works for NULL detection function",{
  mod1_nodf <-dsm(N~s(x,y), NULL, mexdolphins$segdata, mexdolphins$obsdata,
                  strip.width=8000)
  set.seed(1123)
  mod1.movblk_nodf <- dsm.var.movblk(mod1_nodf, mexdolphins$preddata, n.boot=2,
                              block.size=3, off.set=444*1000*1000, bar=FALSE)

  expect_that(summary(mod1.movblk_nodf)$cv[1],
              equals(0.3117824, tol=cv.tol))

  # throw an error if you want detection function uncertainty with no
  # detection function
  expect_error(dsm.var.movblk(mod1_nodf, mexdolphins$preddata, n.boot=2,
                              block.size=3, off.set=444*1000*1000, bar=FALSE,
                              ds.uncertainty=TRUE),
               "Cannot incorporate detection function uncertainty with no detection function!")

})
