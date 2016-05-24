context("Test GetDefaultForm")

test_that("Previous Ys are included in default formulas when they should be", {
  n <- 100
  data <- data.frame(W = rbinom(n, 1, .5),
                     A1 = rbinom(n, 1, .5),
                     Y2 = rbinom(n, 1, .5),
                     A2 = rbinom(n, 1, .5),  
                     Y3 = rbinom(n, 1, .5))

  data2 <- transform(data, Y3 = ifelse(Y2==1, 1, Y3))

  Anodes <- c("A1", "A2") 
  Ynodes <- c("Y2", "Y3")

  forms <- ltmle(data, Anodes, Ynodes=Ynodes, abar=c(1,1), survivalOutcome=FALSE)$formulas
  expect_that("Y2" %in% RhsVars(forms$gform["A2"]), is_true())
  expect_that("Y2" %in% RhsVars(forms$Qform["Y3"]), is_true())  

  forms2 <- ltmle(data2, Anodes, Ynodes=Ynodes, abar=c(1,1), survivalOutcome=TRUE)$formulas
  expect_that("Y2" %in% RhsVars(forms2$gform["A2"]), is_false())
  expect_that("Y2" %in% RhsVars(forms2$Qform["Y3"]), is_false())  
})
