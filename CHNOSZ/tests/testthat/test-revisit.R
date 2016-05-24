context("revisit")

# initial setup
suppressMessages({
  basis("CHNOS")
  basis("O2", -65)
  species(c("leucine", "glycine", "glutamic acid"))
  # a 0-D case
  a <- affinity()
  e0 <- equilibrate(a)
  # a 1-D case
  a <- affinity(O2=c(-75, -65))
  e1 <- equilibrate(a)
  # a 2-D case
  a <- affinity(O2=c(-75, -65), NH3=c(-4, 100, 3))
  e2 <- equilibrate(a)
  # a 3-D case
  a <- affinity(O2=c(-75, -65, 4), H2O=c(-8, 0, 3), NH3=c(-6, -4, 2))
  e3 <- equilibrate(a)
})

test_that("inconsistent arguments produce an error", {
  expect_error(get.objfun("affinity"), "affinity is not a function with an attribute named 'optimum'")
  expect_error(revisit(list(1, 1), plot.it=TRUE), "can't make a plot if 'eout' is not the output from equilibrate\\(\\)")
  expect_error(revisit(list(1, c(1, 2))), "the list provided in 'eout' is not the output from equilibrate\\(\\)")
  expect_error(revisit(e1, "RMSD"), "loga2 must be supplied for RMSD")
  expect_error(revisit(e1, "RMSD", list(1, 1)), "loga2 has different length \\(2\\) than list in eout \\(3\\)")
  # commented because a plot is still initialized ...
  #expect_error(revisit(e2, "CV", style.2D="xxx"), "2D plot style xxx not one of 'contour' or 'image'")
})

test_that("0-D, 1-D, 2-D and 3-D calculations give identical results at the same conditions", {
  r0.qqr <- revisit(e0, "qqr", plot.it=FALSE)
  r1.qqr <- revisit(e1, "qqr", plot.it=FALSE)
  r2.qqr <- revisit(e2, "qqr", plot.it=FALSE)
  r3.qqr <- revisit(e3, "qqr", plot.it=FALSE)
  # check that we get the same values
  expect_equal(c(r0.qqr$H), tail(r1.qqr$H, 1))
  expect_equal(c(r0.qqr$H), r3.qqr$H[4, 3, 2])
  # check that we get the same index and same optimum
  expect_equal(r1.qqr$ixopt, r2.qqr$ixopt, check.attributes=FALSE)
  expect_equal(r1.qqr$optimum, r2.qqr$optimum)
})

test_that("non-referenced objectives give expected results", {
  # the non-referenced objectives use only the logarithms of activities in eout (loga1)
  r1.cv <- revisit(e1, "CV", plot.it=FALSE)
  r1.sd <- revisit(e1, "SD", plot.it=FALSE)
  r1.shannon <- revisit(e1, "shannon", plot.it=FALSE)
  r1.qqr <- revisit(e1, "qqr", plot.it=FALSE)
  # the tests will alert us to significant numerical changes
  # but so far haven't been independently verified
  expect_equal(r1.cv$optimum, 0.30576, tolerance=1e-5) 
  expect_equal(r1.sd$optimum, 0.000284694, tolerance=1e-5) 
  expect_equal(r1.shannon$optimum, 1.066651, tolerance=1e-5)
  expect_equal(r1.qqr$optimum, 0.999783, tolerance=1e-5)
})

test_that("referenced objectives give expected results", {
  # the referenced objectives compare the logarithms of activities (loga1) to reference values (loga2)
  # the spearman correlation coefficient
  r1.spearman <- revisit(e1, "spearman", c(1, 2, 3), plot.it=FALSE)
  expect_equal(head(r1.spearman$H, 1), -1) # perfect anti-rank correlation
  expect_equal(max(r1.spearman$H), 1)  # perfect rank correlation
  # where logarithm of activity of the 3rd species (glutamic acid) maximizes
  r1.logact <- revisit(e1, "logact", 3, plot.it=FALSE)
  expect_equal(r1.logact$ixopt, 71)
})

test_that("DGtr objective gives zero at equilibrium and >0 not at equilibrium", {
  # let's use n-alkanes
  basis(c("CH4", "H2"), c("gas", "gas"))
  species(c("methane", "ethane", "propane", "n-butane"), "liq")
  # calculate equilibrium distribution over a range of logaH2
  a1 <- affinity(H2=c(-10, -5, 101), exceed.Ttr=TRUE)
  e1 <- equilibrate(a1)
  # take the equilibrium distribution at logfH2 = -7.5 as the reference distribution
  loga2 <- list2array(e1$loga.equil)[51, ]
  # calculate the DGtr/RT relative to the reference distribution
  r1 <- revisit(e1, "DGtr", loga2=loga2, plot.it=FALSE)
  # we should find a minimum of zero at logfH2 = -7.5
  expect_equal(min(r1$H), 0)
  expect_equal(r1$xopt, -7.5)
  # we can even go into 2 dimensions
  # (it's a slightly longer test, so don't run it on CRAN)
  if(!any(grepl("R_CHECK_TIMINGS", names(Sys.getenv())))) {
    a2 <- affinity(H2=c(-10, -5, 101), T=c(0, 100, 101), exceed.Ttr=TRUE)
    e2 <- equilibrate(a2)
    r2 <- revisit(e2, "DGtr", loga2=loga2, plot.it=FALSE)
    # we should DGtr=0 at the temperature of the reference distribution (25 degC)
    expect_equal(min(r2$H), 0)
    expect_equal(r2$yopt, 25)
  }
})
