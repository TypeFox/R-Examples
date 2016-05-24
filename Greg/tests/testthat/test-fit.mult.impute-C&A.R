############################
# Tests for the            #
# fit.mult.impute          #
# in printCrudeAndAdjusted #
# issues                   #
############################

context("Tests for fit.mult.impute")

suppressMessages(library("rms"))

test_that("Check regular linear regression with getC&A",{
  set.seed(10)
  n <- 500
  ds <- data.frame(
    ftime = rexp(n),
    fstatus = sample(0:1, size = n, replace = TRUE),
    x1 = factor(sample(LETTERS[1:4], size = n, replace = TRUE)),
    x2 = rnorm(n, mean = 3, 2),
    x3 = factor(sample(letters[1:3], size = n, replace = TRUE)),
    subsetting = factor(sample(c(TRUE, FALSE), size = n, replace = TRUE)))
  
  ds$missing_var_1 <- factor(sample(letters[1:4], size=n, replace=TRUE))
  ds$missing_var_2 <- factor(sample(letters[1:4], size=n, replace=TRUE))
  ds$y <- rnorm(nrow(ds)) +
    (as.numeric(ds$x1)-1) * 1 +
    (as.numeric(ds$missing_var_1)-1)*1 + 
    (as.numeric(ds$missing_var_2)-1)*.5
  
  # Create a messy missing variable
  non_random_missing <- sample(which(ds$missing_var_1 %in% c("b", "d")), 
                               size = 150, replace=FALSE)
  # Restrict the non-random number on the x2 variables
  non_random_missing <- non_random_missing[non_random_missing %in%
                                             which(ds$x2 > mean(ds$x2)*1.5) &
                                             non_random_missing %in%
                                             which(ds$x2 > mean(ds$y))]
  ds$missing_var_1[non_random_missing] <- NA
  
  # Simple missing
  ds$missing_var_2[sample(1:nrow(ds), size=50)] <- NA
  
  # Setup the rms environment
  ddist <<- datadist(ds)
  options(datadist = "ddist")
  
  fit_ols <- ols(y ~ x1 + x2 + x3 + 
                   missing_var_1 + missing_var_2, data=ds) 
  impute_formula <- 
    as.formula(paste("~",
                     paste(all.vars(as.formula(fit_ols)),
                           collapse="+")))
  
  sink(file=ifelse(Sys.info()["sysname"] == "Windows",
                   "NUL",
                   "/dev/null"))
  imp_ds <- aregImpute(impute_formula, data = ds, n.impute = 10)
  
  fmult <- suppressWarnings(fit.mult.impute(formula(fit_ols), 
                                            fitter = lm, xtrans = imp_ds, data = ds))
    
  a <- suppressWarnings(getCrudeAndAdjustedModelData(fmult))
  sink()
  expect_true(sum(a[,"Adjusted"] - coef(fmult)) < .Machine$double.eps)
  
  sink(file=ifelse(Sys.info()["sysname"] == "Windows",
                   "NUL",
                   "/dev/null"))
  fmult <- suppressWarnings(fit.mult.impute(formula(fmult), 
                            fitter = ols, xtrans = imp_ds, data = ds))
  
  a <- suppressWarnings(getCrudeAndAdjustedModelData(fmult))
  sink()
  # Remove the intercept as the fitter was ols
  expect_true(sum(a[,"Adjusted"] - coef(fmult)[-1]) < .Machine$double.eps)
  
  sink(file=ifelse(Sys.info()["sysname"] == "Windows",
                   "NUL",
                   "/dev/null"))
  single_fit <- suppressWarnings(fit.mult.impute(y ~ missing_var_1, 
                                                 fitter = ols, xtrans = imp_ds, data = ds))
  sink()
  expect_true(sum(a[grep("missing_var_1", rownames(a)),"Crude"] - 
                    coef(single_fit)[grep("missing_var_1", 
                                          names(coef(single_fit)))]) < 
                .Machine$double.eps)
  
})

test_that("Check regular linear regression with printC&A",{
  set.seed(10)
  n <- 500
  ds <- data.frame(
    ftime = rexp(n),
    fstatus = sample(0:1, size = n, replace = TRUE),
    x1 = factor(sample(LETTERS[1:4], size = n, replace = TRUE)),
    x2 = rnorm(n, mean = 3, 2),
    x3 = factor(sample(letters[1:3], size = n, replace = TRUE)),
    subsetting = factor(sample(c(TRUE, FALSE), size = n, replace = TRUE)))
  
  ds$missing_var_1 <- factor(sample(letters[1:4], size=n, replace=TRUE))
  ds$missing_var_2 <- factor(sample(letters[1:4], size=n, replace=TRUE))
  ds$y <- rnorm(nrow(ds)) +
    (as.numeric(ds$x1)-1) * 1 +
    (as.numeric(ds$missing_var_1)-1)*1 + 
    (as.numeric(ds$missing_var_2)-1)*.5
  
  # Create a messy missing variable
  non_random_missing <- sample(which(ds$missing_var_1 %in% c("b", "d")), 
                               size = 150, replace=FALSE)
  # Restrict the non-random number on the x2 variables
  non_random_missing <- non_random_missing[non_random_missing %in%
                                             which(ds$x2 > mean(ds$x2)*1.5) &
                                             non_random_missing %in%
                                             which(ds$x2 > mean(ds$y))]
  ds$missing_var_1[non_random_missing] <- NA
  
  # Simple missing
  ds$missing_var_2[sample(1:nrow(ds), size=50)] <- NA
  
  # Setup the rms environment
  ddist <<- datadist(ds)
  options(datadist = "ddist")
  
  fit_ols <- ols(y ~ x1 + x2 + x3 + 
                   missing_var_1 + missing_var_2, data=ds) 
  impute_formula <- 
    as.formula(paste("~",
                     paste(names(attr(fit_ols$terms, "dataClasses")),
                           collapse="+")))
  
  sink(file=ifelse(Sys.info()["sysname"] == "Windows",
                   "NUL",
                   "/dev/null"))
  imp_ds <- aregImpute(impute_formula, data = ds, n.impute = 10)
  
  fmult <- suppressWarnings(fit.mult.impute(formula(fit_ols), 
                                            lm, imp_ds, data = ds))
  
  a <- suppressWarnings(printCrudeAndAdjustedModel(fmult))
  sink()
  expect_equivalent(tail(attr(a, "rgroup"), 1), "missing_var_2")
  
  expect_equivalent(tail(rownames(a), 4),
                    c(tail(levels(ds$missing_var_1), 1), 
                      levels(ds$missing_var_2)[-1]),
                    info=paste("The reference should not be included by default",
                               "without specifying add_reference"))
  
  expect_equivalent(ncol(a), 4)
  
  sink(file=ifelse(Sys.info()["sysname"] == "Windows",
                   "NUL",
                   "/dev/null"))
  a <- suppressWarnings(printCrudeAndAdjustedModel(fmult, add_references = TRUE))
  sink()
  expect_equivalent(tail(rownames(a), 5),
                    c(tail(levels(ds$missing_var_1), 1), 
                      levels(ds$missing_var_2)),
                    info=paste("The reference should be included ",
                               "if add_reference == TRUE"))
  expect_equivalent(ncol(a), 4)
  
  sink(file=ifelse(Sys.info()["sysname"] == "Windows",
                   "NUL",
                   "/dev/null"))
  coef_change_tst <- 
    suppressWarnings(printCrudeAndAdjustedModel(fmult, 
                                                add_references = TRUE, 
                                                impute_args = list(coef_change=TRUE)))
  coef_change_tst_name <- 
    suppressWarnings(printCrudeAndAdjustedModel(fmult, 
                                                add_references = TRUE, 
                                                impute_args = list(coef_change=list(name="cc test"))))
  vi_infl_tst_name <- 
    suppressWarnings(printCrudeAndAdjustedModel(fmult, 
                                                add_references = TRUE, 
                                                impute_args = list(variance.inflation=list(name="vi test"))))
  sink()
  expect_equivalent(ncol(coef_change_tst), 5, 
                    info = "When using fit.mult.impute additional columns should be added")
  
  expect_equivalent(tail(colnames(coef_change_tst_name), 1),
                    "cc test",
                    info = "When using fit.mult.impute additional columns should be nameable")
  
  expect_equivalent(tail(colnames(vi_infl_tst_name), 1),
                    "vi test",
                    info = "When using fit.mult.impute additional columns should be nameable")
  
})

test_that("Check logistic, survival, and poissong for antilog and more",{
  set.seed(10)
  n <- 500
  ds <- data.frame(
    ftime = rexp(n),
    fstatus = sample(0:1, size = n, replace = TRUE),
    x1 = factor(sample(LETTERS[1:4], size = n, replace = TRUE)),
    x2 = rnorm(n, mean = 3, 2),
    x3 = factor(sample(letters[1:3], size = n, replace = TRUE)),
    subsetting = factor(sample(c(TRUE, FALSE), size = n, replace = TRUE)))
  
  ds$missing_var_1 <- factor(sample(letters[1:4], size=n, replace=TRUE))
  ds$missing_var_2 <- factor(sample(letters[1:4], size=n, replace=TRUE))
  ds$y <- rnorm(nrow(ds)) +
    (as.numeric(ds$x1)-1) * 1 +
    (as.numeric(ds$missing_var_1)-1)*1 + 
    (as.numeric(ds$missing_var_2)-1)*.5
  
  # Create a messy missing variable
  non_random_missing <- sample(which(ds$missing_var_1 %in% c("b", "d")), 
                               size = 150, replace=FALSE)
  # Restrict the non-random number on the x2 variables
  non_random_missing <- non_random_missing[non_random_missing %in%
                                             which(ds$x2 > mean(ds$x2)*1.5) &
                                             non_random_missing %in%
                                             which(ds$x2 > mean(ds$y))]
  ds$missing_var_1[non_random_missing] <- NA
  
  # Simple missing
  ds$missing_var_2[sample(1:nrow(ds), size=50)] <- NA
  
  # Setup the rms environment
  ddist <<- datadist(ds)
  options(datadist = "ddist")
  
  s <- Surv(ds$ftime, ds$fstatus == 1)
  sink(file=ifelse(Sys.info()["sysname"] == "Windows",
                   "NUL",
                   "/dev/null"))
  
  imp_ds <- aregImpute(as.formula(paste("~",
                                        paste(colnames(ds),
                                              collapse="+"))),
                       data = ds, n.impute = 10)
  
  fit <- suppressWarnings(fit.mult.impute(s ~ x1 + x2 + strat(x3) + 
                                            missing_var_1 + missing_var_2, 
                                          fitter = cph, data=ds, xtrans = imp_ds))
  
  fit <- coxph(s ~ x1 + x2 + strata(x3) + 
                 missing_var_1 + missing_var_2, data=ds)
  suppressWarnings(printCrudeAndAdjustedModel(fit))
  sink()
  
  fit <- glm(fstatus ~ offset(log(ftime)) + x1 + x2 + x3, data=ds, family = binomial)  
  
  # TODO: add tests
})

