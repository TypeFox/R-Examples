library(testthat)
library(rpf)
#options(error = utils::recover)

context("lmp ICC")

test_that("LMP", {

  ## Check vs. 2PL from same package
  lmp.k0<-rpf.lmp(k=0)
  twopl<-rpf.drm()
  par<-c(.7,-.1)
  par2<-c(exp(.7),-.1,logit(0),logit(1))
  expect_equal(rpf.prob(lmp.k0, par, c(-1,0,1)), rpf.prob(twopl, par2, c(-1,0,1)))

  ## 3rd order polynomial
  lmp.k1 <- rpf.lmp(k=1)

  ## Response probabilities manually copied from code used for Falk & Cai
  lmp.k1.values<-matrix(c(.86737623446,.52497918748,.00005929256,
                          .1326238,.4750208,.99994070744),
                        byrow=TRUE,nrow=2)
  par<-c(.7,-.1,-2,.5)

  ## Since not obtaining above from external code, there will be some minor
  ## differences if we compare to too many decimal places
  expect_equal(lmp.k1.values,rpf.prob(lmp.k1, par, c(-1,0,1)), tolerance=1e-7)

  ## 5th order polynomial
  lmp.k2 <- rpf.lmp(k=2)

  ## Response probabilities manually copied from code used for Falk & Cai
  lmp.k2.values<-matrix(c(.9998936,.6224593,.4339060,
                          .0001063994,.3775406677,.5660940386),
                        byrow=TRUE,nrow=2)
  par<-c(.6,-.5,.4,-.2,1,-1)

  ## Since not obtaining above from external code, there will be some minor
  ## differences if we compare to too many decimal places
  expect_equal(lmp.k2.values,rpf.prob(lmp.k2, par, c(-1,0,1)), tolerance=1e-7)

})
