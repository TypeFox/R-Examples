library(testthat)

set.seed(10)
n <- 200
cov <- data.frame(
  ftime = rexp(n),
  fstatus = sample(0:1, n, replace=TRUE),
  x1 = runif(200),
  x2 = runif(200),
  x3 = runif(200),
  x_factor = sample(LETTERS[1:4], n, replace = TRUE))

library(rms)
ddist <<- datadist(cov)
options(datadist = "ddist")

context("Check fortestplotCombineRegrObj")
test_that("Check regular lm", {
  #TODO: Add tests - current code is only for coverage reasons

  fit1 <- cph(Surv(ftime, fstatus) ~ x1 + x2, data=cov)
  fit2 <- cph(Surv(ftime, fstatus) ~ x1 + x3, data=cov)
  
  forestplotCombineRegrObj (
    regr.obj = list(fit1, fit2),
    variablesOfInterest.regexp = "(x2|x3)",
    reference.names = c("First model", "Second model"),
    new_page = TRUE)
  
  modifyNameFunction <- function(x){
    if (x == "x1")
      return ("Covariate A")
    
    if (x == "x2")
      return (expression(paste("My ", beta[2])))
    
    return (x)
  }
  
  forestplotCombineRegrObj (
    regr.obj = list(fit1, fit2),
    variablesOfInterest.regexp = "(x2|x3)",
    reference.names = c("First model", "Second model"),
    rowname.fn = modifyNameFunction,
    new_page = TRUE)
})

test_that("Test getModelData4Forestplot",{
  # simulated data to test 
  fit1 <- cph(Surv(ftime, fstatus) ~ x1 + x2, data=cov)
  fit2 <- cph(Surv(ftime, fstatus) ~ x1 + x3, data=cov)
  
  # regr.obj, 
  # exp = TRUE, 
  # variablesOfInterest.regexp,
  # ref_labels,
  # add_first_as_ref
  data <- getModelData4Forestplot(
    regr.obj = list(fit1, fit2),
    exp = TRUE,
    variablesOfInterest.regexp = "(x2|x3)",
    ref_labels = FALSE,
    add_first_as_ref = FALSE)
  expect_equivalent(data[[1]]["x2","beta"],
                    exp(coef(fit1)["x2"]))
  expect_equivalent(data[[2]]["x3","beta"],
                    exp(coef(fit2)["x3"]))
  expect_false("x1" %in% rownames(data[[1]]))
  expect_false("x1" %in% rownames(data[[2]]))

  data <- getModelData4Forestplot(
    regr.obj = list(fit1, fit2),
    exp = TRUE,
    variablesOfInterest.regexp = "(x1)",
    ref_labels = FALSE,
    add_first_as_ref = FALSE)
  expect_equivalent(data[[1]]["x1","beta"],
                    exp(coef(fit1)["x1"]))
  expect_equivalent(data[[2]]["x1","beta"],
                    exp(coef(fit2)["x1"]))
  expect_true("x1" %in% rownames(data[[1]]))
  expect_true("x1" %in% rownames(data[[2]]))

  
  fit3 <- cph(Surv(ftime, fstatus) ~ x1 + x2 + x_factor, data=cov)
  fit4 <- cph(Surv(ftime, fstatus) ~ x1 + x3 + x_factor, data=cov)
  
  data <- getModelData4Forestplot(
    regr.obj = list(fit3, fit4),
    exp = TRUE,
    variablesOfInterest.regexp = "(x1|x_factor)",
    ref_labels = FALSE,
    add_first_as_ref = FALSE)
  expect_equivalent(data[[1]]["x1","beta"],
                    exp(coef(fit3)["x1"]))
  expect_equivalent(data[[2]]["x1","beta"],
                    exp(coef(fit4)["x1"]))
  expect_equivalent(data[[1]]["x_factor=B","beta"],
                    exp(coef(fit3)["x_factor=B"]))
  expect_equivalent(data[[2]]["x_factor=C","beta"],
                    exp(coef(fit4)["x_factor=C"]))
  expect_true("x1" %in% rownames(data[[1]]))
  expect_true("x1" %in% rownames(data[[2]]))
  expect_false("x2" %in% rownames(data[[1]]))
  expect_false("x2" %in% rownames(data[[2]]))
  expect_equivalent(nrow(data[[1]]),
                    4)
  expect_equivalent(nrow(data[[2]]),
                    4)
  
  data <- getModelData4Forestplot(
    regr.obj = list(fit3, fit4),
    exp = TRUE,
    variablesOfInterest.regexp = "(x_factor)",
    ref_labels = "x_factor=A",
    add_first_as_ref = TRUE)
  expect_equivalent(data[[1]][1,"beta"],
                    1)
  expect_equivalent(data[[2]][1,"beta"],
                    1)

  data <- getModelData4Forestplot(
    regr.obj = list(m1 = fit3, m2 = fit4),
    exp = FALSE,
    variablesOfInterest.regexp = "(x_factor)",
    ref_labels = "x_factor=A",
    add_first_as_ref = TRUE)
  expect_equivalent(data[[1]][1,"beta"],
                    0)
  expect_equivalent(data[[2]][1,"beta"],
                    0)
  
  expect_equivalent(data[[1]]["x_factor=B","beta"],
                    coef(fit3)["x_factor=B"])
  expect_equivalent(data[[2]]["x_factor=C","beta"],
                    coef(fit4)["x_factor=C"])
})
