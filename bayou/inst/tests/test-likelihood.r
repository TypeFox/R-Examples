context("can calculate likelihoods")
test_that("can calculate likelihoods", {
  data(chelonia)
  tree <- chelonia$phy
  dat <- chelonia$dat
  cache <- bayou:::.prepare.ou.univariate(tree, dat, SE=0)
  pars <- list(alpha=0.01, sig2=1, k=3, ntheta=4, theta=c(3,4,5,6), sb= c(408, 399, 448), loc=c(8, 9, 31), t2=2:4)
  bayou.lik(QGpars, cache, dat, model="QG")$loglik
  expect_that(is.finite(OU.lik(pars, tree, dat)$loglik), is_true())
  expect_that(is.finite(bayou.lik(pars, cache, dat)$loglik), is_true())
  expect_that(bayou.lik(pars, cache, dat)$loglik[1], equals(OU.lik(pars, tree, dat)$loglik[1]))
  expect_that(bayou.lik(pars, cache, dat)$loglik[1], equals(-494.775911175547))
  ##Test Brownian motion works
  pars <- list(alpha=0, sig2=1, k=3, ntheta=4, theta=c(3,4,5,6), sb= c(408, 399, 448), loc=c(8, 9, 31), t2=2:4)
  geiger_bm <- bm.lik(tree, dat, model="BM")
  expect_that(OU.lik(pars, cache, dat)$loglik[1], equals(geiger_bm(c(pars$sig2, 0, pars$theta[1]), root=ROOT.GIVEN)[1]))
  expect_that(bayou.lik(pars, cache, dat)$loglik[1], equals(geiger_bm(c(pars$sig2, 0, pars$theta[1]), root=ROOT.GIVEN)[1]))
  ##Test QG parameterization
  QGpars <- list(h2=0.4, P=2, Ne=100, w2=10, k=3, ntheta=4, theta=c(3,4,5,6), sb= c(408, 399, 448), loc=c(8, 9, 31), t2=2:4)
  expect_that(bayou.lik(QGpars, cache, dat, model="QG")$loglik[1], equals(OU.lik(QGpars, tree, dat, model="QG")$loglik[1]))
  ##Test OUrepar parameterization
  OUrepars <- list(halflife=50, Vy=3, k=3, ntheta=4, theta=c(3,4,5,6), sb= c(408, 399, 448), loc=c(8, 9, 31), t2=2:4)
  expect_that(bayou.lik(OUrepars, cache, dat, model="OUrepar")$loglik[1], equals(OU.lik(OUrepars, tree, dat, model="OUrepar")$loglik[1]))
  ##Deprecated tests
  #expect_that(is.finite(emOU.lik(pars,emap,cache,dat,SE=0.1,model="OU")$loglik[1]),is_true())
  #expect_that(emOU.lik(pars,emap,tree,dat,SE=0.1,model="OU")$loglik[1],
  #            equals(emOU.lik(pars,emap,cache,dat,SE=0.1,model="OU",method="invert")$loglik))
  #expect_that(emOU.lik(pars,emap,cache,dat,SE=0.1,model="OU")$loglik[1],
  #            equals(emOU.lik(pars,emap,cache,dat,SE=0.1,model="OU",method="invert")$loglik))
  #expect_that(.emOU.lik(pars,emap,cache,dat,SE=0.1,model="OU")$loglik[1],
  #            equals(emOU.lik(pars,emap,cache,dat,SE=0.1,model="OU",method="invert")$loglik))
  #expect_that(emOU.lik(QGpars,emap,cache,dat,SE=0.1,model="QG")$loglik[1],
  #            equals(emOU.lik(pars,emap,cache,dat,SE=0.1,model="OU",method="invert")$loglik))
  #expect_that(.emOU.lik(QGpars,emap,cache,dat,SE=0.1,model="QG")$loglik[1],
  #            equals(emOU.lik(pars,emap,cache,dat,SE=0.1,model="OU",method="invert")$loglik))
  #expect_that(smOU.lik(pars,tree,dat,SE=0.1,model="OU")$loglik[1],
  #            equals(emOU.lik(pars,emap,cache,dat,SE=0.1,model="OU",method="invert")$loglik))
  #expect_that(smOU.lik(QGpars,tree,dat,SE=0.1,model="QG")$loglik[1],
  #            equals(emOU.lik(pars,emap,cache,dat,SE=0.1,model="OU",method="invert")$loglik))
  #expect_that(.smOU.lik(pars,cache,dat,SE=0.1,model="OU")$loglik[1],
  #            equals(emOU.lik(pars,emap,cache,dat,SE=0.1,model="OU",method="invert")$loglik))
  #expect_that(.smOU.lik(QGpars,cache,dat,SE=0.1,model="QG")$loglik[1],
  #            equals(emOU.lik(pars,emap,cache,dat,SE=0.1,model="OU",method="invert")$loglik))
  #pars$sb <- which(emap$sh==1)
  #pars$loc <- emap$r1[emap$sh==1]
  #pars$t2 <- emap$t2[emap$sh==1]
  #QGpars$sb <- which(emap$sh==1)
  #QGpars$loc <- emap$r1[emap$sh==1]
  #QGpars$t2 <- emap$t2[emap$sh==1]
  #expect_that(.OU.lik(pars,cache,dat,SE=0.1,model="OU")$loglik[1],
  #            equals(emOU.lik(pars,emap,cache,dat,SE=0.1,model="OU",method="invert")$loglik))
  #expect_that(OU.lik(pars,tree,dat,SE=0.1,model="OU")$loglik[1],
  #            equals(emOU.lik(pars,emap,cache,dat,SE=0.1,model="OU",method="invert")$loglik))
  #expect_that(.OU.lik(QGpars,cache,dat,SE=0.1,model="QG")$loglik[1],
  #            equals(emOU.lik(pars,emap,cache,dat,SE=0.1,model="OU",method="invert")$loglik))
  })
