context("Excursions.inla")

test_that("excursions.inla, test ind", {
  if(require.nowarnings("INLA")){
  data <- testdata.inla()
  ind1 = c(1,2,3,4)
  ind2 = c(4,3,2,1)
  ind3 = rep(FALSE,data$n)
  ind3[1:4] = TRUE

  res1 = excursions.inla(data$result, data$stack, ind=ind1,method="QC",
                          tag="pred", u=0, type='>', seed = data$seed)
  res2 = excursions.inla(data$result, data$stack, ind=ind2,method="QC",
                          tag="pred", u=0, type='>', seed = data$seed)
  res3 = excursions.inla(data$result, data$stack, ind=ind3,method="QC",
                          tag="pred", u=0, type='>', seed = data$seed)

  expect_equal(res1$F,res2$F,tolerance=1e-7)
  expect_equal(res2$F,res3$F,tolerance=1e-7)
  }
})
