context("Contourmap.inla")

test_that("Contourmap.inla, test ind", {
  if (require.nowarnings("INLA")) {
  data <- testdata.inla()
  ind1 = c(1,2,3,4)
  ind2 = c(4,3,2,1)
  ind3 = rep(FALSE,data$n)
  ind3[1:4] = TRUE

  res1 = contourmap.inla(data$result, data$stack, tag = "pred",
                         n.levels=2,ind=ind1, seed=data$seed,alpha=0.1)
  res2 = contourmap.inla(data$result, data$stack, tag = "pred",
                         n.levels=2,ind=ind2, seed=data$seed,alpha=0.1)
  res3 = contourmap.inla(data$result, data$stack, tag = "pred",
                         n.levels=2,ind=ind3, seed=data$seed,alpha=0.1)

  expect_equal(res1$F,res2$F,tolerance=1e-7)
  expect_equal(res2$F,res3$F,tolerance=1e-7)
}
})
