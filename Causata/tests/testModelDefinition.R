# Test of the ModelDefinition class
# Author: David Barker <support@causata.com>

library(testthat)
library(glmnet)
library(Causata)

context("ModelDefinition")

equals <- testthat::equals

# This is just a simple test that the refactoring into ModelDefinition and VariableDefintion works.
#
test_that("PMML Generation from glmnet model with no transformations", {
  iv <-  1:20
  dv <- 21:40 + rnorm(5)
  df <- data.frame(iv__AP=iv)
  data <- CausataData(df, dv)
  
  formula.obj = formula(dependent.variable ~ iv__AP)
  model <- cv.glmnet(family="gaussian", x=model.matrix(formula.obj, data$df), y=dv, alpha=0.9, nfolds=3)
  model.def <- ModelDefinition(model, data, formula.obj)
  
  variable.def <- VariableDefinition("my-model")
  pmml <- ToPmml(model.def, variable.def)
  
  expect_that(is.null(pmml), equals(FALSE))
})
