
context("solarPar")

### basic examples

test_that("kurtosis", {
  data(dat30)
  mod <- solarPolygenic(trait1 ~ 1, dat30)

  val <- solarParKurtosis(mod)
  expect_true(!is.na(val))
})

test_that("covlist pvalue", {
  data(dat30)

  mod0 <- solarPolygenic(trait1 ~ sex, dat30)
  mod1 <- solarPolygenic(trait1 ~ sex, dat30, covtest = TRUE)
  mod2 <- solarPolygenic(trait1 ~ age + sex, dat30, covtest = TRUE)

  val0 <- solarParCovlistP(mod0)
  val1 <- solarParCovlistP(mod1)
  val2 <- solarParCovlistP(mod2)

  expect_true(is.na(val0)) # p-value is not computed if covtest is FALSE
  expect_true(is.numeric(val1)) 
  expect_equal(length(val1), 1) # a single value is expected, as one covariate is specified in `mod1`
  expect_true(is.numeric(val2))
  expect_equal(length(val2), 2) # two values are expected, as two covariates are specified in `mod2`
})


### derived functions

test_that("function explainedVarProp", {
  data(dat30)
  mod0 <- solarPolygenic(trait1 ~ 1, dat30)
  tab0 <- explainedVarProp(mod0)
  
  mod1 <- solarPolygenic(trait1 ~ age + age^2 + sex, dat30)
  tab1 <- explainedVarProp(mod1)

  mod2 <- solarPolygenic(trait1 ~ age + age^2 + sex, dat30, covtest = TRUE)
  tab2 <- explainedVarProp(mod2)
  
  expect_true(all(is.na(tab0$explainedVarProp))) # no results for a model without any covariate
  
  expect_true(any(is.na(tab1$explainedVarProp))) # no matter whether covtest is TRUE/FALSE, the results must be for the case of all covariates 
  expect_equal(sum(!is.na(tab1$explainedVarProp)), 1) 
  
  expect_true(all(!is.na(tab2$explainedVarProp))) # if covtest is TRUE, all results must be obtained
})

