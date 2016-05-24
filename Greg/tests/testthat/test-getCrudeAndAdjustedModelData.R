library("testthat")
library("survival")

set.seed(10)
n <- 500
ds <<- data.frame(
  ftime = rexp(n),
  fstatus = sample(0:1, size = n, replace = TRUE),
  y = rnorm(n = n),
  x1 = factor(sample(LETTERS[1:4], size = n, replace = TRUE)),
  x2 = rnorm(n, mean = 3, 2),
  x3 = factor(sample(letters[1:3], size = n, replace = TRUE)),
  boolean = sample(0L:1L, size = n, replace = TRUE),
  subsetting = factor(sample(c(TRUE, FALSE), size = n, replace = TRUE)))

library(rms)
ddist <<- datadist(ds)
options(datadist="ddist")

context("getCrudeAndAdjustedModelData")
test_that("Correct number of rows and columns", {
  fit1 <- coxph(Surv(ds$ftime, ds$fstatus == 1) ~ x1 + x2 + x3, data=ds)
  
  # Check that it doesn't include the rcs() spline since this doesn't
  # make sense 
  library(splines)
  fit2 <- coxph(Surv(ds$ftime, ds$fstatus == 1) ~ x1 + ns(x2, 4) + strata(x3), data=ds)
  
  data_matrix <- getCrudeAndAdjustedModelData(fit1)
  expect_equivalent(NROW(data_matrix), 3 + 1 + 2)
  expect_equivalent(NCOL(data_matrix), 6)
  
  data_matrix <- getCrudeAndAdjustedModelData(fit2)
  expect_equivalent(NROW(data_matrix), 3 + 0 + 0)
  expect_equivalent(NCOL(data_matrix), 6)
})

test_that("Same order of rows and matching results", {
  fit1 <- coxph(Surv(ds$ftime, ds$fstatus == 1) ~ x1 + x2 + x3, data=ds)
  data_matrix <- getCrudeAndAdjustedModelData(fit1)
  expect_equal(rownames(data_matrix), names(coef(fit1)))
  expect_true(sum(data_matrix[,"Adjusted"] - exp(coef(fit1))) <
                .Machine$double.eps)
  expect_true(sum(data_matrix[,tail(1:ncol(data_matrix), 2)] -
                    exp(confint(fit1))) < .Machine$double.eps)
  
  unadjusted_fit <- coxph(Surv(ds$ftime, ds$fstatus == 1) ~ x2, data=ds)
  expect_true(data_matrix["x2","Crude"] == exp(coef(unadjusted_fit)))
  expect_true(all(data_matrix["x2",2:3] == exp(confint(unadjusted_fit))))
})

test_that("Check subsetting", {
  fit1 <- coxph(Surv(ds$ftime, ds$fstatus == 1) ~ x1 + x2 + x3, data=ds, 
                subset=subsetting == TRUE)
  data_matrix <- getCrudeAndAdjustedModelData(fit1)
  expect_equal(rownames(data_matrix), names(coef(fit1)))
  expect_true(sum(data_matrix[,"Adjusted"] -
                    exp(coef(fit1))) < .Machine$double.eps)
  expect_true(sum(data_matrix[,tail(1:ncol(data_matrix), 2)] - 
                    exp(confint(fit1))) < .Machine$double.eps)
  
  unadjusted_fit <- update(fit1, .~x2)
  expect_true(data_matrix["x2","Crude"] == exp(coef(unadjusted_fit)))
  expect_true(all(data_matrix["x2",2:3] == exp(confint(unadjusted_fit))))

  # Doesn't seem to work to get the with() environment
#   fit1 <- with(ds, coxph(Surv(ftime, fstatus == 1) ~ x1 + x2 + x3,  
#                          subset=subsetting == TRUE))
#   data_matrix <- getCrudeAndAdjustedModelData(fit1)
#   expect_equal(rownames(data_matrix), names(coef(fit1)))
#   expect_equal(data_matrix[,"Adjusted"], exp(coef(fit1)))
#   expect_equal(data_matrix[,tail(1:ncol(data_matrix), 2)], exp(confint(fit1)))

  attach(ds)
  fit1 <- coxph(Surv(ftime, fstatus == 1) ~ x1 + x2 + x3,  
                         subset=subsetting == TRUE)
  data_matrix <- getCrudeAndAdjustedModelData(fit1)
  expect_equal(rownames(data_matrix), names(coef(fit1)))
  expect_true(sum(data_matrix[,"Adjusted"] - 
                    exp(coef(fit1))) < .Machine$double.eps)
  expect_true(sum(data_matrix[,tail(1:ncol(data_matrix), 2)] -
                    exp(confint(fit1))) < .Machine$double.eps)

  unadjusted_fit <- update(fit1, .~x2)
  expect_true(data_matrix["x2","Crude"] == exp(coef(unadjusted_fit)))
  expect_true(all(data_matrix["x2",2:3] == exp(confint(unadjusted_fit))))

  unadjusted_fit <- update(fit1, .~x2, subset =  subsetting == FALSE)
  expect_false(data_matrix["x2","Crude"] == exp(coef(unadjusted_fit)))
  detach(ds)
})

test_that("Same order of rows", {
  fit1 <- cph(Surv(ds$ftime, ds$fstatus == 1) ~ x1 + x2 + x3, data=ds)
  data_matrix <- getCrudeAndAdjustedModelData(fit1)
  expect_match(rownames(data_matrix)[1:(length(levels(ds$x1))-1)], "^x1")
  expect_match(rownames(data_matrix)[length(levels(ds$x1))], "^x2")
  expect_match(rownames(data_matrix)[(length(levels(ds$x1)) +1):nrow(data_matrix)], "^x3 - [bc]")
  
  unadjusted_fit <- cph(Surv(ds$ftime, ds$fstatus == 1) ~ x2, data=ds)
  expect_true(data_matrix["x2","Crude"] - exp(coef(unadjusted_fit)) < .Machine$double.eps)
  expect_true(all(data_matrix["x2",2:3] - exp(confint(unadjusted_fit)) < .Machine$double.eps))
  
})

test_that("A few bug tests",{
  
  # Produced an  error with integer as input due to sort:
  #   Error in value.chk(at, i, factors[[jf]], 0, Limval) : 
  #     character value not allowed for variable boolean
  fit1 <- cph(Surv(ds$ftime, ds$fstatus == 1) ~ x1 + x2 + x3 + boolean, data=ds)
  expect_is(getCrudeAndAdjustedModelData(fit1), "matrix")
  
})

test_that("Same order of rows", {
  fit_lm <- lm(y ~ x1 + x2 + x3, data=ds) 
  
  data_matrix <- getCrudeAndAdjustedModelData(fit_lm)
  expect_match(rownames(data_matrix)[1], "^\\(Intercept\\)")
  expect_match(rownames(data_matrix)[2], "^x1")
  expect_match(rownames(data_matrix)[5], "^x2")
  expect_match(rownames(data_matrix)[6:7], "^x3[bc]")

})

test_that("Correct values for rows - lm", {
  fit_lm <- lm(y ~ x1 + x2 + x3, data=ds) 
  
  data_matrix <- getCrudeAndAdjustedModelData(fit_lm)
  expect_true(all(data_matrix[,"Adjusted"] ==coef(fit_lm)))
  expect_true(all(data_matrix[,tail(1:ncol(data_matrix), 2)] ==confint(fit_lm)))
  
  fit_lm_ua <- lm(y ~ x3, data=ds)
  expect_true(all(data_matrix[grep("^x3", rownames(data_matrix)),
                              "Crude"] ==
                    coef(fit_lm_ua)[2:3]))
  expect_true(all(data_matrix[grep("^x3", rownames(data_matrix)),
                              2:3] ==
                    confint(fit_lm_ua)[2:3,]))
  
  fit_lm_ua <- lm(y ~ x1, data=ds)
  expect_true(all(data_matrix[grep("^x1", rownames(data_matrix)),
                              "Crude"] ==
                    coef(fit_lm_ua)[-1]))
  expect_true(all(data_matrix[grep("^x1", rownames(data_matrix)),
                              2:3] ==
                    confint(fit_lm_ua)[-1,]))
})

test_that("Correct values for rows - ols", {
  fit_ols <- ols(y ~ x1 + x2 + x3, data=ds) 
  
  # The rms-package does not report the intercept in summary
  # and we therefore skip
  data_matrix <- getCrudeAndAdjustedModelData(fit_ols)
  expect_true(sum(data_matrix[,"Adjusted"] -
                    coef(fit_ols)[-1]) < .Machine$double.eps)
  expect_true(sum(data_matrix[,tail(1:ncol(data_matrix), 2)] -
                    confint(fit_ols)[-1,])  < .Machine$double.eps)
  
  fit_ols_ua <- ols(y ~ x3, data=ds)
  data_matrix <- getCrudeAndAdjustedModelData(fit_ols_ua)
  expect_true(sum(data_matrix[grep("^x3", rownames(data_matrix)), "Crude"] - 
                    coef(fit_ols_ua)[2:3])  < .Machine$double.eps)
  expect_true(sum(data_matrix[grep("^x3", rownames(data_matrix)), 2:3] -
                    confint(fit_ols_ua)[2:3,]) < .Machine$double.eps)
  
  fit_ols_ua <- ols(y ~ x1, data=ds)
  data_matrix <- getCrudeAndAdjustedModelData(fit_ols_ua)
  expect_true(sum(data_matrix[grep("^x1", rownames(data_matrix)), "Crude"] - 
                    coef(fit_ols_ua)[-1]) < .Machine$double.eps)
  expect_true(sum(data_matrix[grep("^x1", rownames(data_matrix)), 2:3] - 
                    confint(fit_ols_ua)[-1,]) < .Machine$double.eps)
})

test_that("Correct values for rows - glm", {
  fit_glm <- glm(y ~ x1 + x2 + x3, data=ds) 
  
  data_matrix <- getCrudeAndAdjustedModelData(fit_glm)
  expect_true(all(abs(data_matrix[,"Adjusted"] - 
                        coef(fit_glm)) < .Machine$double.eps))
  expect_true(all(abs(data_matrix[,tail(1:ncol(data_matrix), 2)] -
                        confint(fit_glm)) < .Machine$double.eps))
  
  fit_glm_ua <- glm(y ~ x3, data=ds)
  expect_true(all(abs(data_matrix[grep("^x3", rownames(data_matrix)), "Crude"] -
                        coef(fit_glm_ua)[2:3]) < .Machine$double.eps))
  expect_true(all(abs(data_matrix[grep("^x3", rownames(data_matrix)), 2:3] -
                        confint(fit_glm_ua)[2:3,]) < .Machine$double.eps ))
  
  fit_glm_ua <- glm(y ~ x1, data=ds)
  expect_true(all(abs(data_matrix[grep("^x1", rownames(data_matrix)), "Crude"] -
                        coef(fit_glm_ua)[-1]) < .Machine$double.eps))
  expect_true(all(abs(data_matrix[grep("^x1", rownames(data_matrix)), 2:3] -
                    confint(fit_glm_ua)[-1,]) < .Machine$double.eps))
})

test_that("Variable selection - default", {
  fit_glm <- glm(y ~ x1 + x2 + x3, data=ds) 
  
  data_matrix <- getCrudeAndAdjustedModelData(fit_glm, var_select = c("x1"))
  expect_true(sum(data_matrix[,"Adjusted"] -
                    coef(fit_glm)[grep("^x1", names(coef(fit_glm)))]) < .Machine$double.eps)
  
  fit_glm_ua <- glm(y ~ x1, data=ds)
  expect_true(sum(data_matrix[, "Crude"] - 
                    coef(fit_glm_ua)[-1]) < .Machine$double.eps)
  
  data_matrix <- getCrudeAndAdjustedModelData(fit_glm, var_select = c("Intercept", "x1"))
  expect_true(sum(data_matrix[,"Adjusted"] -
                    coef(fit_glm)[grep("(^x1|ntercept)", names(coef(fit_glm)))]) < .Machine$double.eps)
  
  fit_glm_ua <- glm(y ~ x1, data=ds)
  expect_true(sum(data_matrix[-1, "Crude"] -
                    coef(fit_glm_ua)[-1]) < .Machine$double.eps)

  expect_error(getCrudeAndAdjustedModelData(fit_glm, var_select = c("x222222")))
})

test_that("Variable selection - rms", {
  fit_glm <- Glm(y ~ x1 + x2 + x3, data=ds) 
  
  data_matrix <- getCrudeAndAdjustedModelData(fit_glm, var_select = c("x1"))
  expect_true(sum(data_matrix[,"Adjusted"] - 
                    coef(fit_glm)[grep("^x1", names(coef(fit_glm)))]) < .Machine$double.eps)
  
  fit_glm_ua <- update(fit_glm, . ~ x1)
  expect_true(sum(data_matrix[, "Crude"] - 
                    coef(fit_glm_ua)[-1]) < .Machine$double.eps)
  
  data_matrix <- getCrudeAndAdjustedModelData(fit_glm, var_select = c("x2", "x1"))
  expect_true(sum(data_matrix[,"Adjusted"] -
                    coef(fit_glm)[grep("(^x[12])", names(coef(fit_glm)))]) < .Machine$double.eps)
  
  expect_error(getCrudeAndAdjustedModelData(fit_glm, var_select = c("x222222")))
})

test_that("Strata and clusters should be still present in the crude format coxph",{
  test1 <- data.frame(time=c(4,3,1,1,2,2,3), 
                      status=c(1,1,1,0,1,1,0), 
                      x1=c(0,2,1,1,1,0,0), 
                      x2=c(1,3,5,2,0,9,1), # Additional to the coxph example 
                      sex=c(0,0,0,0,1,1,1)) 

  # Create the simplest test data set 
  # Fit a stratified model 
  fit <- coxph(Surv(time, status) ~ x1 + x2 + strata(sex), data = test1) 
  x <- getCrudeAndAdjustedModelData(fit)
  
  fit <- update(fit, . ~ x1 + strata(sex)) 
  expect_true(sum(x["x1", "Crude"] - 
                    exp(coef(fit))) < .Machine$double.eps)

  fit <- update(fit, . ~ x2 + strata(sex)) 
  expect_true(sum(x["x2", "Crude"] - 
                    exp(coef(fit))) < .Machine$double.eps)
  
  fit <- coxph(Surv(time, status) ~ x1 + x2 + strata(sex), data = test1) 
  x <- getCrudeAndAdjustedModelData(fit, remove_strata = TRUE)
  fit <- update(fit, . ~ x1) 
  expect_true(sum(x["x1", "Crude"] -
                    exp(coef(fit))) < .Machine$double.eps)
  
  fit <- update(fit, . ~ x2) 
  expect_true(sum(x["x2", "Crude"] - 
                    exp(coef(fit))) < .Machine$double.eps)
  
  fit <- coxph(Surv(time, status) ~ x1 + x2 + cluster(sex), data = test1) 
  x <- getCrudeAndAdjustedModelData(fit, remove_cluster = FALSE)
  fit <- update(fit, . ~ x1 + cluster(sex)) 
  expect_true(all(x["x1", 2:3] -exp(confint(fit)) < .Machine$double.eps))
  
  fit <- update(fit, . ~ x2 + cluster(sex)) 
  expect_true(all(x["x2", 2:3] -exp(confint(fit)) < .Machine$double.eps))
  
  fit <- coxph(Surv(time, status) ~ x1 + x2 + cluster(sex), data = test1) 
  x <- getCrudeAndAdjustedModelData(fit, remove_cluster = TRUE)
  fit <- update(fit, . ~ x1) 
  expect_true(all(x["x1", 2:3] - exp(confint(fit)) < .Machine$double.eps))
  
  fit <- update(fit, . ~ x2) 
  expect_true(all(x["x2", 2:3] -exp(confint(fit)) < .Machine$double.eps))
})

test_that("Strata and clusters should be still present in the crude format cph",{
  library(rms)
  test1 <<- data.frame(time=c(4,3,1,1,2,2,3), 
                       status=c(1,1,1,0,1,1,0), 
                       x1=c(0,2,1,1,1,0,0), 
                       x2=c(1,3,5,2,0,9,1), # Additional to the coxph example 
                       sex=c(0,0,0,0,1,1,1)) 
  ddist <<- datadist(test1)
  options(datadist="ddist")
  fit <- cph(Surv(time, status) ~ x1 + x2 + strat(sex), 
             data = test1) 
  x <- getCrudeAndAdjustedModelData(fit)
  
  fit <- update(fit, . ~ x1 + strat(sex)) 
  expect_true(sum(x["x1", "Crude"] - 
                    exp(coef(fit))) < .Machine$double.eps)
  
  fit <- update(fit, . ~ x2 + strat(sex)) 
  expect_true(sum(x["x2", "Crude"] - 
                    exp(coef(fit))) < .Machine$double.eps)
  
  fit <- cph(Surv(time, status) ~ x1 + x2 + strat(sex), data = test1) 
  x <- getCrudeAndAdjustedModelData(fit, remove_strata = TRUE)
  fit <- update(fit, . ~ x1) 
  expect_true(sum(x["x1", "Crude"] - 
                    exp(coef(fit))) < .Machine$double.eps)
  
  fit <- update(fit, . ~ x2) 
  expect_true(sum(x["x2", "Crude"] - 
                    exp(coef(fit))) < .Machine$double.eps)
  
  #########################
  # Same but for clusters #
  # - note that clusters  #
  #   only affect confid. #
  #   intervals           #
  #########################
  
  fit <- cph(Surv(time, status) ~ x1 + x2 + cluster(sex), 
             data = test1)
  x <- getCrudeAndAdjustedModelData(fit, remove_cluster = FALSE)
  fit <- update(fit, . ~ x1 + cluster(sex)) 
  expect_true(all(x["x1", 2:3] -exp(confint(fit)) < .Machine$double.eps))
  
  fit <- update(fit, . ~ x2 + cluster(sex)) 
  expect_true(all(x["x2", 2:3] -exp(confint(fit)) < .Machine$double.eps))
})