library(dsm)
library(Distance)
library(testthat)

par.tol<-1e-5

context("test inputs")

test_that("formula specs",{

  # load the Gulf of Mexico dolphin data
  data(mexdolphins)

  # fit a detection function and look at the summary
  hn.model <- ds(mexdolphins$distdata, max(mexdolphins$distdata$distance),
                 adjustment = NULL)

  ## models for count
  count.gcv <- 42.9169051
  count.N<-dsm(N~s(x,y), hn.model, mexdolphins$segdata, mexdolphins$obsdata)
  expect_that(count.N$gcv.ubre, equals(count.gcv,tolerance=par.tol))

  count.n<-dsm(n~s(x,y), hn.model, mexdolphins$segdata, mexdolphins$obsdata)
  expect_that(count.n$gcv.ubre, equals(count.gcv,tolerance=par.tol))

  count.count<-dsm(count~s(x,y), hn.model, mexdolphins$segdata,
                   mexdolphins$obsdata)
  expect_that(count.count$gcv.ubre, equals(count.gcv,tolerance=par.tol))

  count.abundance<-dsm(abundance~s(x,y), hn.model, mexdolphins$segdata,
                       mexdolphins$obsdata)
  expect_that(count.abundance$gcv.ubre, equals(count.gcv,tolerance=par.tol))


  ## models for abund.est
  abund.est.gcv <- 57.3159048
  abund.est.Nhat<-dsm(Nhat~s(x,y), hn.model, mexdolphins$segdata,
                      mexdolphins$obsdata)
  expect_that(abund.est.Nhat$gcv.ubre, equals(abund.est.gcv,tolerance=par.tol))
  abund.est.abund.est<-dsm(abundance.est~s(x,y), hn.model, mexdolphins$segdata,
                           mexdolphins$obsdata)
  expect_that(abund.est.abund.est$gcv.ubre,equals(abund.est.gcv,
                                                  tolerance=par.tol))
  abund.est.abundance<-dsm(abundance.est~s(x,y), hn.model, mexdolphins$segdata,
                       mexdolphins$obsdata)
  expect_that(abund.est.abundance$gcv.ubre, equals(abund.est.gcv,
                                                   tolerance=par.tol))


  ## models for density
  D.gcv <- 1.660703e-07
  density.D<-dsm(D~s(x,y), hn.model, mexdolphins$segdata, mexdolphins$obsdata,
                 weights=rep(1,nrow(mexdolphins$segdata)))
  expect_that(density.D$gcv.ubre, equals(D.gcv,tolerance=par.tol))

  density.density<-dsm(density~s(x,y), hn.model, mexdolphins$segdata,
                       mexdolphins$obsdata,
                       weights=rep(1,nrow(mexdolphins$segdata)))
  expect_that(density.density$gcv.ubre, equals(D.gcv,tolerance=par.tol))

  density.Dhat<-dsm(Dhat~s(x,y), hn.model, mexdolphins$segdata,
                    mexdolphins$obsdata,
                    weights=rep(1,nrow(mexdolphins$segdata)))
  expect_that(density.Dhat$gcv.ubre, equals(D.gcv,tolerance=par.tol))

  density.density.est<-dsm(density.est~s(x,y), hn.model, mexdolphins$segdata,
                           mexdolphins$obsdata,
                           weights=rep(1,nrow(mexdolphins$segdata)))
  expect_that(density.density.est$gcv.ubre, equals(D.gcv,tolerance=par.tol))

  # check that Effort is not zero
  mex_zero_effort <- mexdolphins$segdata
  mex_zero_effort$Effort[c(1,5,10)] <- 0
  expect_error(dsm(abundance.est~s(x,y), hn.model, mex_zero_effort,
                           mexdolphins$obsdata))

})

test_that("Missing columns cause errors",{

  data(mexdolphins)

  seg <- mexdolphins$segdata
  obs <- mexdolphins$obsdata

  for(mcov in c("object","Sample.Label","size","distance")){
    obs_missing <- mexdolphins$obsdata
    obs_missing[[mcov]] <- NULL
    expect_error(dsm(N~s(x,y), NULL, seg, obs_missing, segment.area = 8000^2),
                 paste0("Column(s) \"",mcov,
                        "\" not found in observation.data."),
                 fixed=TRUE)
  }

  for(mcov in c("Effort","Sample.Label")){
    seg_missing <- seg
    seg_missing[[mcov]] <- NULL
    expect_error(dsm(N~s(x,y), NULL, seg_missing, obs, segment.area = 8000^2),
                 paste0("Column(s) \"",mcov,
                        "\" not found in segment.data."),
                 fixed=TRUE)
  }


})



