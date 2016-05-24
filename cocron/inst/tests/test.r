library(testthat)

context("cocron.n.coefficients")

test_that("Calculate examples from the literature", {
  # Feldt, L. S., Woodruff, D. J., Salih, F. A. (1987). Statistical inference for coefficient alpha. Applied Psychological Measurement, 11, 93-103.

  result <- cocron.n.coefficients(alpha=c(.784,.875,.936), items=c(5,5,5), n=c(51,101,151), dep=FALSE)
  expect_that(result@df, is_identical_to(2))
  expect_that(round(result@statistic,3), is_identical_to(22.926))


  r <- rbind(
  c(1,.8,.6,.75),
  c(NA,1,.65,.7),
  c(NA,NA,1,.55),
  c(NA,NA,NA,1)
  )
  result <- cocron.n.coefficients(alpha=c(.857,.875,.800,.833), items=c(50,40,35,25), n=100, dep=TRUE, r=r)
  expect_that(result@df, is_identical_to(3))
  expect_that(round(result@statistic,2), is_identical_to(10.66))
  expect_that(round(result@p.value,3), is_identical_to(.014))
})

test_that("Compare with results from AlphaTest", {
  # Lautenschlager, G. J., & Meade, A. W., 2008. AlphaTest

  #1
  result <- cocron.n.coefficients(alpha=c(.6,.7), items=c(50,100), n=c(150,200), dep=FALSE)
  expect_that(result@df, is_identical_to(1))
  expect_that(round(result@statistic,4), is_identical_to(3.4535))
  expect_that(round(result@p.value,4), is_identical_to(.0631))

  #2
  r <- rbind(
  c(1,.8),
  c(NA,1)
  )
  result <- cocron.n.coefficients(alpha=c(.6,.5), items=c(35,30), n=90, dep=TRUE, r=r)
  expect_that(result@df, is_identical_to(1))
  expect_that(round(result@statistic,4), is_identical_to(2.8699))
  expect_that(round(result@p.value,4), is_identical_to(.0903))

  #3
  result <- cocron.n.coefficients(alpha=c(.9,.7), items=c(25,30), n=c(200,150), dep=FALSE)
  expect_that(result@df, is_identical_to(1))
  expect_that(round(result@statistic,3), is_identical_to(47.559))
  expect_that(round(result@p.value,4), is_identical_to(.0000))

  #4
  r <- rbind(
  c(1,.5),
  c(NA,1)
  )
  result <- cocron.n.coefficients(alpha=c(.5,.4), items=c(20,20), n=c(200, 100), dep=TRUE, r=r)
  expect_that(result@df, is_identical_to(1))
  expect_that(round(result@statistic,4), is_identical_to(1.4895))
  expect_that(round(result@p.value,4), is_identical_to(.2223))

  #5
  result <- cocron.n.coefficients(alpha=c(.9,.7,.8), items=c(15,15,12), n=c(100,150,130), dep=FALSE)
  expect_that(result@df, is_identical_to(2))
  expect_that(round(result@statistic,4), is_identical_to(31.3795))
  expect_that(round(result@p.value,4), is_identical_to(.0000))

  #6
  result <- cocron.n.coefficients(alpha=c(.7,.75,.7), items=c(10,10,10), n=c(100,120,130), dep=FALSE)
  expect_that(result@df, is_identical_to(2))
  expect_that(round(result@statistic,4), is_identical_to(1.0553))
  expect_that(round(result@p.value,4), is_identical_to(.59))

  #7
  r <- rbind(
  c(1,.5,.2),
  c(NA,1,.8),
  c(NA,NA,1)
  )
  result <- cocron.n.coefficients(alpha=c(.75,.7,.72), items=c(12,12,12), n=c(60,70,80), dep=TRUE, r=r)
  expect_that(result@df, is_identical_to(2))
  expect_that(round(result@statistic,4), is_identical_to(0.7130))
  expect_that(round(result@p.value,4), is_identical_to(.7001))

  #8
  r <- rbind(
  c(1,.3,.2),
  c(NA,1,.7),
  c(NA,NA,1)
  )
  result <- cocron.n.coefficients(alpha=c(.7,.75,.72), items=c(13,14,15), n=c(140,130,110), dep=TRUE, r=r)
  expect_that(result@df, is_identical_to(2))
  expect_that(round(result@statistic,4), is_identical_to(1.1747)) # the result from AlphaTest is actually 1.1732
  expect_that(round(result@p.value,4), is_identical_to( 0.5558)) # the result from AlphaTest is actually .5562

  #9
  result <- cocron.n.coefficients(alpha=c(.7,.8,.75,.77), items=c(15,15,20,20), n=c(100,100,200,200), dep=FALSE)
  expect_that(result@df, is_identical_to(3))
  expect_that(round(result@statistic,4), is_identical_to(4.9738))
  expect_that(round(result@p.value,4), is_identical_to(.1737))

  #10
  r <- rbind(
  c(1,.5,.6,.7),
  c(NA,1,.4,.2),
  c(NA,NA,1,.2),
  c(NA,NA,NA,1)
  )
  result <- cocron.n.coefficients(alpha=c(.75,.72,.70,.77), items=c(15,15,15,15), n=c(100,100,150,150), dep=TRUE, r=r)
  expect_that(result@df, is_identical_to(3))
  expect_that(round(result@statistic,4), is_identical_to(2.9236))
  expect_that(round(result@p.value,4), is_identical_to(.4036))
})

context("cocron.two.coefficients")

test_that("Calculate examples from the literature", {
  # Charter, R. A., Feldt, L. S. (1996). Testing the equality of two alpha coefficients. Perceptual and Motor Skills, 82, 763-768.

  result <- cocron.two.coefficients(alpha=c(.78,.71), n=c(41,151), dep=FALSE, alternative="greater")
  expect_that(result@df, is_identical_to(c(40, 150)))
  expect_that(round(result@statistic,3), is_identical_to(1.318))
  expect_that(round(result@p.value,3) * 2, is_identical_to(.242)) # did not use alternative="two.sided" because it seems that the p-value reported in the paper has been rounded before it was doubled

  result <- cocron.two.coefficients(alpha=c(.82,.89), n=27, dep=TRUE, r=.74, alternative="two.sided")
  expect_that(result@df, is_identical_to(25))
  expect_that(round(result@statistic,3), is_identical_to(1.849))
})


context("cronbach.alpha.CI")

test_that("Calculate examples from the literature", {
  # Feldt, L. S., Woodruff, D. J., Salih, F. A. (1987). Statistical inference for coefficient alpha. Applied Psychological Measurement, 11, 93-103. See p. 95 formulas 8 and 9.

  expect_that(round(cronbach.alpha.CI(alpha=.790, n=41, items=26, conf.level=.9),3), equals(c(lower.bound=.705,upper.bound=.862))) # The actual values reported in the article are L = .704 and U = .861, however, the authors used F values that are only accurate to two decimal places.
})
