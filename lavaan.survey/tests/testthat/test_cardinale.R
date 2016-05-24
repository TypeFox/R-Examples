
library(lavaan.survey)

context("Cardinale example")

data(cardinale)

fit.card <- sem('
   PatchDiversity ~ logNutrient + logNutrient2 + StreamDiversity
   Biomass ~ PatchDiversity + logNutrient
   O2Production ~ logNutrient + Biomass
   logNutrient ~~ logNutrient2', data = cardinale, 
                fixed.x = FALSE, estimator = "MLM")


des.card <- svydesign(ids = ~Stream, probs = ~1, data = cardinale)
fit.card.survey <- lavaan.survey(fit.card, des.card, estimator = "MLM")
fit.card.survey.MLMVS <- lavaan.survey(fit.card, des.card, estimator = "MLMVS")

test_that("scaled chi square matches (MLM)", {
  fm <- fitMeasures(fit.card.survey)
  
  expect_equal(fm['chisq.scaling.factor'], 1.519, tolerance = 1e-3, check.attributes = FALSE)
  expect_equal(fm['chisq'], 5.901, tolerance = 1e-3, check.attributes = FALSE)
  expect_equal(fm['df'], 7, check.attributes = FALSE)
  
})

test_that("scaled chi sq matches (MLMVS)", {
  fm <- fitMeasures(fit.card.survey.MLMVS)
  
  expect_equal(fm['chisq.scaling.factor'], 4.593529, tolerance = 1e-4, check.attributes = FALSE)
  expect_equal(fm['df.scaled'], 2.316, tolerance = 1e-3, check.attributes = FALSE)
  expect_equal(fm['pvalue.scaled'], 0.600, tolerance = 1e-3, check.attributes = FALSE)
})

test_that("F test matches", {
  skip_on_cran() # May not have CompQuadForm
  
  expect_equal(pval.pFsum(fit.card.survey, survey.design = des.card), 0.609748, tolerance = 1e-5, check.attributes = FALSE)
  expect_equal(pval.pFsum(fit.card.survey, survey.design = des.card, method="integration"), 0.6231761, tolerance = 1e-5, check.attributes = FALSE)
  expect_equal(pval.pFsum(fit.card.survey, survey.design = des.card, method="satterthwaite")[1], 0.6079676, tolerance = 1e-5, check.attributes = FALSE)
})

test_that("standard errors match", {
  ses <- sqrt(diag(vcov(fit.card.survey)))
  expect_equal(ses, c(0.39374, 0.16843, 0.08605, 0.00734, 0.03419, 4e-04, 0.00425, 0.20201, 18.04856, 0.06281, 8e-05, 0.0741, 0.30242, 189.73026, 6.84561, 0.24455, 0.00505, 0.03056, 0.07347, 5.6909), tolerance = 1e-5, check.attributes = FALSE)
  
})

test_that("an estimate matches", {
 expect_equal(coef(fit.card.survey)[1], 0.3619211, tolerance = 1e-5, check.attributes = FALSE) 
})