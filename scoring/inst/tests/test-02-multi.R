context('multicategory')

test_that('manyalt', {
  data(WeatherProbs)
  scores <- calcscore(tcat ~ tblw + tnrm + tabv, fam="pow",
                      param=c(2,rep(1/3,3)), data=WeatherProbs,
                      bounds=c(0,1))
  expect_is(scores, 'numeric')
  expect_false(any(is.na(scores)))
  expect_true(all(scores <= 1) & all(scores >= 0))

  scores2 <- calcscore(tcat ~ tblw + tnrm + tabv, fam="sph",
                       param=c(2,rep(1/3,3)), data=WeatherProbs,
                       bounds=c(0,1), reverse=TRUE)
  expect_is(scores2, 'numeric')
  expect_false(any(is.na(scores2)))
  expect_true(all(scores2 <= 1) & all(scores2 >= 0))

  scores3 <- calcscore(tcat ~ tblw + tnrm + tabv, fam="sph",
                       param=c(3,.1,.1,.8), data=WeatherProbs,
                       bounds=c(0,1), reverse=TRUE)
  expect_is(scores3, 'numeric')
  expect_false(any(is.na(scores3)))
  expect_true(all(scores3 <= 1) & all(scores3 >= 0))

  scores4 <- calcscore(tcat ~ tblw + tnrm + tabv, fam="sph",
                       param=c(3,.1,.1,.8), data=WeatherProbs[1:500,],
                       bounds=c(0,1), reverse=TRUE, ordered=TRUE)
  expect_is(scores4, 'numeric')
  expect_false(any(is.na(scores4)))
  expect_true(all(scores4 <= 1) & all(scores4 >= 0))

  ## Fixed bug in ordered scores, when outcome is in the
  ## "final" column of forecasts.  This check ensures that
  ## the better forecast receives the higher score (because
  ## reverse=TRUE).
  expect_true(scores4[53] > scores4[60])

  r2 <- seq(0, .6, .05)
  r <- cbind(.4, r2, .6 - r2)
  j <- rep(1, length(r2))
  quad <- calcscore(j ~ r, fam="pow", param=2, bounds=c(-1,1), reverse=TRUE)
  expect_true(all((quad >= -1) & (quad <= 1)))

})
