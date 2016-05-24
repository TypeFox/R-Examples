library("testthat")
context("isFit")

library("survival")
library("rms")

test_that("Cox regression isFit works", {
  set.seed(10)
  n <- 500
  ds <<- data.frame(
    ftime = rexp(n),
    fstatus = sample(0:1, size = n, replace = TRUE),
    y = rnorm(n = n),
    x1 = factor(sample(LETTERS[1:4], size = n, replace = TRUE)),
    x2 = rnorm(n, mean = 3, 2),
    x3 = factor(sample(letters[1:3], size = n, replace = TRUE)))
  
  dd <<- datadist(ds)
  options(datadist = "dd")
  
  s <- Surv(ds$ftime, ds$fstatus == 1)
  fit1 <- coxph(s ~ x1 + x2 + x3, data=ds)
  
  s <- Surv(ds$ftime, ds$fstatus == 1)
  fit_cox <- coxph(s ~ x1 + x2 + x3, data=ds)
  fit_cph <- cph(s ~ x1 + x2 + x3, data=ds)
  
  fit_logistic <- glm(fstatus ~ x1 + x2 + x3, data=ds, family=binomial)
  fit_lrm <- lrm(fstatus ~ x1 + x2 + x3, data=ds)

  expect_true(isFitCoxPH(fit_cox))
  
  expect_true(isFitCoxPH(fit_cph))
  
  expect_false(isFitCoxPH(fit_logistic))

  expect_false(isFitCoxPH(fit_lrm))
})


test_that("Logit regression isFit works", {
  set.seed(10)
  n <- 500
  ds <- data.frame(
    ftime = rexp(n),
    fstatus = sample(0:1, size = n, replace = TRUE),
    y = rnorm(n = n),
    x1 = factor(sample(LETTERS[1:4], size = n, replace = TRUE)),
    x2 = rnorm(n, mean = 3, 2),
    x3 = factor(sample(letters[1:3], size = n, replace = TRUE)))
  
  dd <<- datadist(ds)
  options(datadist = "dd")
  
  s <- Surv(ds$ftime, ds$fstatus == 1)
  fit1 <- coxph(s ~ x1 + x2 + x3, data=ds)
  
  s <- Surv(ds$ftime, ds$fstatus == 1)
  fit_cox <- coxph(s ~ x1 + x2 + x3, data=ds)
  fit_cph <- cph(s ~ x1 + x2 + x3, data=ds)
  
  fit_logistic <- glm(fstatus ~ x1 + x2 + x3, data=ds, family=binomial)
  fit_lrm <- lrm(fstatus ~ x1 + x2 + x3, data=ds)
  
  expect_false(isFitLogit(fit_cox))
  
  expect_false(isFitLogit(fit_cph))
  
  expect_true(isFitLogit(fit_logistic))
  
  expect_true(isFitLogit(fit_lrm))
})
