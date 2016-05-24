context("weight matrix can be calculated")
test_that("weight matrix can be calculated", {
  data(chelonia.simmap)
  tree <- chelonia.simmap$tree
  dat <- chelonia.simmap$dat
  emap <- chelonia.simmap$emap
  cache <- .prepare.ou.univariate(tree, dat)
  pars <- list(alpha=0.1, sig2=1, k=16, optima=c(3,4,5,6), ntheta=4)
  TotExp <- exp(-cache$height*pars$alpha)
  expect_that(apply(simmap.W(tree, pars),1,sum),equals(rep(1,226)))
  expect_that(apply(simmap.W(cache, pars),1,sum),equals(rep(1,226)))
  expect_that(apply(.simmap.W(cache,pars),1,sum),equals(rep(1,226)))
  expect_that(apply(edgemap.W(tree,pars,emap),1,sum),equals(rep(1,226)))
  expect_that(apply(.edgemap.W(cache,pars,emap,TotExp),1,sum),equals(rep(1,226)))
  expect_that(simmap.W(tree,pars),equals(edgemap.W(tree,pars,emap)))
  pars$sb <- which(emap$sh==1)
  pars$loc <- emap$r1[emap$sh==1]
  pars$t2 <- emap$t2[emap$sh==1]
  expect_that(.parmap.W(cache,pars), equals(edgemap.W(tree,pars,emap)))
  expect_that(parmap.W(tree,pars), equals(edgemap.W(tree,pars,emap)))
})
