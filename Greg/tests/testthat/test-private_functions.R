library(testthat)

context("prGetModelVariables")
test_that("Check how well variables are identified and stratification, etc are removed", {
  set.seed(10)
  n <- 50
  ds <- data.frame(
    ftime = rexp(n),
    fstatus = sample(0:1, size = n, replace = TRUE),
    a = factor(sample(c("a", "aa", "aaa"), size = n, replace = TRUE)),
    aa = rnorm(n, mean = 3, 2),
    aaa = factor(sample(c("a", "aa", "aaa"), size = n, replace = TRUE)))

  fit <- lm(ftime ~ ., data=ds)
  expect_equivalent(prGetModelVariables(fit),
                    colnames(ds)[-1])
  
  fit <- lm(ftime ~ a, data=ds)
  expect_equivalent(prGetModelVariables(fit),
                    "a")
  
  fit <- lm(ftime ~ a + aa, data=ds)
  expect_equivalent(prGetModelVariables(fit),
                    c("a", "aa"))
  
  fit <- lm(ftime ~ a * aa, data=ds)
  expect_equivalent(prGetModelVariables(fit),
                    c("a", "aa"))
  
  fit <- lm(ftime ~ a * aa, data=ds)
  expect_equivalent(length(prGetModelVariables(fit, 
                                               remove_interaction_vars = TRUE)),
                    0)

  fit <- lm(ftime ~ a * aa + aaa, data=ds)
  expect_equivalent(length(prGetModelVariables(fit, 
                                               remove_interaction_vars = TRUE)),
                    1)

  fit <- lm(ftime ~ a : aa, data=ds)
  expect_equivalent(length(prGetModelVariables(fit)),
                    0)
  
  fit <- lm(ftime ~ a : aa + a, data=ds)
  expect_equivalent(prGetModelVariables(fit),
                    c("a"))
  
  fit <- lm(ftime ~ offset(aa) + a, data=ds)
  expect_equivalent(prGetModelVariables(fit),
                    c("a"))
  
  fit <- glm(ftime ~ I(aa*2) + a, data=ds)
  expect_equivalent(prGetModelVariables(fit),
                    c("a"))
  
  fit <- with(ds, glm(ftime ~ I(aa*2) + a + aa:a))
  expect_equivalent(prGetModelVariables(fit),
                    c("a"))
  
  library(rms)
  ddist <<- datadist(ds)
  options(datadist = "ddist")
  fit <- Glm(fstatus ~ a * aa, data=ds, family=binomial)
  expect_equivalent(prGetModelVariables(fit),
                    c("a", "aa"))
  expect_equivalent(prGetModelVariables(fit, 
                                        remove_interaction_vars = TRUE),
                    character(0))
  
  library(nlme)
  fit <- lme( distance ~ age , data=Orthodont)
  expect_equivalent(prGetModelVariables(fit),
                    c("age"))
  fit <- lme( distance ~ age +Sex, data=Orthodont, random=~1)
  expect_equivalent(prGetModelVariables(fit),
                    c("age", "Sex"))
  
  fit <- nlme(height ~ SSasymp(age, Asym, R0, lrc),
              data = Loblolly,
              fixed = Asym + R0 + lrc ~ 1,
              random = Asym ~ 1,
              start = c(Asym = 103, R0 = -8.5, lrc = -3.3))
  expect_equivalent(prGetModelVariables(fit),
                    c("Asym","R0","lrc"))
  
  # Bug:
#   fit <- ols(ftime ~ ., data=ds)
#   expect_equivalent(prGetModelVariables(fit),
#                     colnames(ds)[-1])
  
  fit <- ols(ftime ~ a, data=ds)
  expect_equivalent(prGetModelVariables(fit),
                    "a")
  
  fit <- ols(ftime ~ a + aa, data=ds)
  expect_equivalent(prGetModelVariables(fit),
                    c("a", "aa"))
  
  fit <- ols(ftime ~ a * aa, data=ds)
  expect_equivalent(prGetModelVariables(fit),
                    c("a", "aa"))
  
  fit <- ols(ftime ~ a * aa + aaa, data=ds)
  expect_equivalent(length(prGetModelVariables(fit, 
                                               remove_interaction_vars = TRUE)),
                    1)

#Bug:
#   fit <- ols(ftime ~ a : aa, data=ds)
#   expect_equivalent(length(prGetModelVariables(fit)),
#                     0)
  
#   fit <- ols(ftime ~ a : aa + a, data=ds)
#   expect_equivalent(prGetModelVariables(fit),
#                     c("a"))
  
  fit <- ols(ftime ~ offset(aa) + a, data=ds)
  expect_equivalent(prGetModelVariables(fit),
                    c("a"))
  
  fit <- Glm(ftime ~ I(aa*2) + a, data=ds)
  expect_equivalent(prGetModelVariables(fit),
                    c("a"))

#Bug:  
#   fit <- Glm(ftime ~ I(aa*2) + a + aa:a, data=ds)
#   expect_equivalent(prGetModelVariables(fit),
#                     c("a"))
  

  fit <- cph(Surv(ftime, fstatus) ~ a + aa + aaa, data=ds)
  expect_equivalent(prGetModelVariables(fit),
                    c("a", "aa", "aaa"))

  fit <- cph(Surv(ftime, fstatus) ~ a + aa + strat(aaa), data=ds)
  expect_equivalent(prGetModelVariables(fit),
                    c("a", "aa"))

  fit <- cph(Surv(ftime, fstatus) ~ a + I(aa^2) + strat(aaa), data=ds)
  expect_equivalent(prGetModelVariables(fit),
                    c("a"))
  
  fit <- coxph(Surv(ftime, fstatus) ~ a + I(aa^2) + strata(aaa), data=ds)
  expect_equivalent(prGetModelVariables(fit),
                    c("a"))
})

context("prMapVariable2Name")
test_that("Check how the mapping of rows work", {
  set.seed(10)
  n <- 500
  ds <- data.frame(
    y = rnorm(n),
    a = factor(sample(c("a", "aa", "aaa"), size = n, replace = TRUE)),
    aa = rnorm(n, mean = 3, 2),
    aaa = factor(sample(c("a", "aa", "aaa"), size = n, replace = TRUE)))
  
  fit <- lm(y~., data=ds)
  expect_error(prMapVariable2Name(prGetModelVariables(fit),
                                  names(coef(fit)),
                                  data=ds),
               info="Impossible to resolve the true relationship")
  
  ds <- data.frame(
    y = rnorm(n),
    a = factor(sample(c("a", "aa", "aaa"), size = n, replace = TRUE)),
    aa = rnorm(n, mean = 3, 2),
    baa = factor(sample(c("a", "aa", "aaa"), size = n, replace = TRUE)))
  fit <- lm(y~., data=ds)
  out <- prMapVariable2Name(prGetModelVariables(fit),
                            names(coef(fit)),
                            data=ds)
  
  expect_equivalent(out$a$location, 2:3)
  expect_equivalent(out$aa$location, 4)
  expect_equivalent(out$baa$location, 5:6)
  
  ds <- data.frame(
    y = rnorm(n),
    a = factor(sample(c("a", "aa", "aaa"), size = n, replace = TRUE)),
    aaa = rnorm(n, mean = 3, 2),
    baa = factor(sample(c("a", "aa", "aaa"), size = n, replace = TRUE)))
  fit <- lm(y~., data=ds)

  expect_error(prMapVariable2Name(prGetModelVariables(fit),
                                  names(coef(fit)),
                                  data=ds))

  ds <- data.frame(
    y = rnorm(n),
    a = factor(sample(c("a", "aa", "aaa"), size = n, replace = TRUE)),
    ada = rnorm(n, mean = 3, 2),
    aal = factor(sample(c("a", "aa", "aaa"), size = n, replace = TRUE)))
  fit <- lm(y~a*ada+aal, data=ds)
  prMapVariable2Name(prGetModelVariables(fit, remove_interaction_vars = TRUE),
                     names(coef(fit)),
                     data=ds)
  
  fit <- lm(y~a+ada+aal, data=ds)
  out <- prMapVariable2Name(prGetModelVariables(fit, 
                                                add_intercept = TRUE,
                                                remove_interaction_vars = TRUE),
                     names(coef(fit)),
                     data=ds)
  expect_equal(out$`(Intercept)`$location, 1)
  expect_equal(out$`a`$location, 2:3)
  expect_equal(out$`ada`$location, 4)
  expect_equal(out$`aal`$location, 5:6)
  
  fit <- ols(y~a+ada+aal, data=ds)
  out <- prMapVariable2Name(prGetModelVariables(fit, 
                                                add_intercept = TRUE,
                                                remove_interaction_vars = TRUE),
                            names(coef(fit)),
                            data=ds)
  expect_equal(out$`Intercept`$location, 1)
  expect_equal(out$`a`$location, 2:3)
  expect_equal(out$`ada`$location, 4)
  expect_equal(out$`aal`$location, 5:6)
  
  ds <- data.frame(
    ftime = rexp(n),
    fstatus = sample(0:1, size = n, replace = TRUE),
    a = factor(sample(c("a", "aa", "aaa"), size = n, replace = TRUE)),
    aa = rnorm(n, mean = 3, 2),
    aaa = factor(sample(c("a", "aa", "aaa"), size = n, replace = TRUE)))
  
  fit <- cph(Surv(ftime, fstatus)~a+aa+aaa, data=ds)
  out  <- prMapVariable2Name(prGetModelVariables(fit, 
                                                 add_intercept = TRUE,
                                                 remove_interaction_vars = TRUE),
                             names(coef(fit)),
                             data=ds)
  expect_equivalent(length(out), 3)
  expect_equal(out$`a`$location, 1:2)
  expect_equal(out$`aa`$location, 3)
  expect_equal(out$`aaa`$location, 4:5)
  
  ds <- data.frame(
    x1 = factor(sample(c("a", "aa", "aaa"), size = n, replace = TRUE)),
    x2 = rnorm(n, mean = 3, 2))
  expect_error(prMapVariable2Name(c("x1", "x2"),
                                  available_names = "x2",
                                  data=ds))
  
  out <- prMapVariable2Name(c("x1", "x2"),
                            available_names = "x2",
                            data=ds, 
                            force_match = FALSE)
  expect_equivalent(length(out),
                    1)
  expect_equivalent(out$x2$location,
                    1)
  expect_equivalent(out$x2$no_rows,
                    1)
  
  ds <- structure(list(y = c(1L, 0L, 0L, 1L, 0L, 0L, 0L, 0L, 1L, 0L), 
                       x1 = structure(c(4L, 4L, 4L, 2L, 1L, 3L, 1L, 1L, 1L, 4L), 
                                      .Label = c("A", 
                                                 "B", "C", "D"), class = "factor"), 
                       x2 = structure(c(1L, 1L, 
                                        2L, 1L, 1L, 2L, 2L, 2L, 1L, 1L), 
                                      .Label = c("No", "Yes"), class = "factor"), 
                       x3 = structure(c(2L, 1L, 1L, 2L, 1L, 2L, 2L, 1L, 2L, 2L), 
                                      .Label = c("No", "Yes"), 
                                      class = "factor")), 
                  .Names = c("y", "x1", "x2", "x3"),
                  row.names = c(NA, 10L), class = "data.frame")
  out <- prMapVariable2Name(var_names = c('(Intercept)', 'x1', 'x2', 'x3'),
                            available_names = c("x1B", "x1C", "x1D", "x2Yes"),
                            data=ds,
                            force_match = FALSE)
  expect_equivalent(length(out),
                    2)
  expect_equivalent(out$x1$location,
                    1:3)
  expect_equivalent(out$x2$location,
                    4)
})


context("prExtractOutcomeFromModel")
test_that("Handling cox regression survival object", {
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
  
  fit <- coxph(Surv(ds$ftime, ds$fstatus == 1) ~ x1 + x2 + x3, data=ds)
  
  expect_equivalent(prExtractOutcomeFromModel(fit), 
                    as.numeric(ds$fstatus == 1))
  
  s <- Surv(ds$ftime, ds$fstatus == 0)
  fit <- coxph(s ~ x1 + x2 + x3, data=ds)
  
  expect_equivalent(prExtractOutcomeFromModel(fit), 
                    as.numeric(ds$fstatus == 0))

  library(rms)
  ddist <<- datadist(ds)
  options(datadist="ddist")
  
  fit <- cph(Surv(ds$ftime, ds$fstatus == 1) ~ x1 + x2 + x3, data=ds)
  
  expect_equivalent(prExtractOutcomeFromModel(fit), 
                    as.numeric(ds$fstatus == 1))
  
  s <- Surv(ds$ftime, ds$fstatus == 0)
  fit <- cph(s ~ x1 + x2 + x3, data=ds)
  
  expect_equivalent(prExtractOutcomeFromModel(fit), 
                    as.numeric(ds$fstatus == 0))
  
})

test_that("Handling simple linear regression outcomes", {
  set.seed(10)
  n <- 500
  ds <- data.frame(
    y = rnorm(n = n),
    y_fact = sample(c("A (1)", 
                      "A (2)",
                      "B"),
                    size = n,
                    replace=TRUE),
    x1 = factor(sample(LETTERS[1:4], size = n, replace = TRUE)),
    x2 = rnorm(n, mean = 3, 2),
    x3 = factor(sample(letters[1:3], size = n, replace = TRUE)),
    boolean = sample(0L:1L, size = n, replace = TRUE),
    subsetting = factor(sample(c(TRUE, FALSE), size = n, replace = TRUE)))

  fit <- lm(y ~ x1 + x2 + x3, data=ds)
  
  expect_equivalent(prExtractOutcomeFromModel(fit), 
                    ds$y)
  
  fit <- lm(ds$y ~ x1 + x2 + x3, data=ds)
  
  expect_equivalent(prExtractOutcomeFromModel(fit), 
                    ds$y)
  

  fit <- glm(ds$y_fact == "A (1)" ~ x1 + x2 + x3, data=ds)
  
  expect_equivalent(prExtractOutcomeFromModel(fit), 
                    ds$y_fact == "A (1)")
  

  fit <- glm(y_fact == "B" ~ x1 + x2 + x3, data=ds)
  expect_equivalent(prExtractOutcomeFromModel(fit), 
                    ds$y_fact == "B")
})

test_that("Subsetting and missing data", {
  set.seed(10)
  n <- 500
  ds <- data.frame(
    y = rnorm(n = n),
    x1 = factor(sample(LETTERS[1:4], size = n, replace = TRUE)),
    x2 = rnorm(n, mean = 3, 2),
    x3 = factor(sample(letters[1:3], size = n, replace = TRUE)),
    subsetting = factor(sample(c(TRUE, FALSE), size = n, replace = TRUE)))
  
  fit <- lm(y ~ x1 + x2 + x3, data=ds, subset = subsetting == TRUE)
  
  expect_equivalent(prExtractOutcomeFromModel(fit), 
                    subset(ds, subsetting == TRUE)$y)

  # Add missing
  ds$x1[sample(1:nrow(ds), size = 50)] <- NA
  fit <- lm(y ~ x1 + x2 + x3, data=ds)
  
  expect_equivalent(prExtractOutcomeFromModel(fit, mf = ds), 
                    subset(ds)$y)
  
  expect_equivalent(length(prExtractOutcomeFromModel(fit)), 
                    nrow(ds) - 50)

  # Advanced test with missing and varible not in the original dataset
  # in the subsetting argument
  local_var <- TRUE
  fit <- lm(y > 0 ~ x1 + x2 + x3, data=ds, subset = subsetting == local_var)
  expect_equivalent(prExtractOutcomeFromModel(fit, mf = ds), 
                    subset(ds)$y > 0)
})


context("prGetModelData tests")
test_that("Subsetting", {
  set.seed(10)
  n <- 10
  ds <- data.frame(
    y = rnorm(n = n),
    x1 = factor(sample(LETTERS[1:4], size = n, replace = TRUE)),
    x2 = rnorm(n, mean = 3, 2),
    x3 = factor(sample(letters[1:3], size = n, replace = TRUE)),
    subsetting = factor(sample(c(TRUE, FALSE), size = n, replace = TRUE)))
  
  fit <- lm(y ~ x1 + x2 + x3, data=ds, subset = subsetting == TRUE)
  for (n in c("y", "x1", "x2", "x3")){
    expect_equivalent(prGetModelData(fit)[,n,drop=FALSE], 
                      subset(ds, subsetting == TRUE, n))
  }
  
  fit <- lm(y ~ ., data=ds)
  fit <- update(fit, .~.-subsetting, subset=subsetting==TRUE)
  for (n in c("y", "x1", "x2", "x3")){
    expect_equivalent(prGetModelData(fit)[,n,drop=FALSE], 
                      subset(ds, subsetting == TRUE, n))
  }
  
  
  fit <- with(ds, lm(y ~ x1 + x2 + x3, subset = subsetting == TRUE))
  for (n in c("y", "x1", "x2", "x3")){
    expect_equivalent(prGetModelData(fit)[,n,drop=FALSE], 
                      subset(ds, subsetting == TRUE, n))
  }

  attach(ds)
  fit <- lm(y ~ x1 + I(x2^2) + x3, subset = subsetting == TRUE)
  for (n in c("y", "x1", "x2", "x3")){
    expect_equivalent(prGetModelData(fit)[,n,drop=FALSE], 
                      subset(ds, subsetting == TRUE, n))
  }
  detach(ds)
  
})
