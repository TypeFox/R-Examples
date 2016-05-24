context("printCrudeAndAdjustedModel")
test_that("Check position of reference", {
  set.seed(10)
  n <- 500
  ds <- data.frame(
    ftime = rexp(n),
    fstatus = sample(0:1, size = n, replace = TRUE),
    x = factor(sample(LETTERS[1:4], size = n, replace = TRUE)),
    same_label = factor(sample(c("Yes", "No"), size = n, replace = TRUE)),
    same_labell = factor(sample(c("Yes", "No"), size = n, replace = TRUE)),
    same_labelll = factor(sample(c("Yes", "No"), size = n, replace = TRUE)),
    boolean = sample(c(TRUE, FALSE), size = n, replace = TRUE),
    subsetting = factor(sample(c(TRUE, FALSE), size = n, replace = TRUE)))
  
  library(survival)
  fit <- coxph(Surv(ftime, fstatus == 1) ~ x + boolean, data=ds)
  
  a <- printCrudeAndAdjustedModel(fit, add_references=TRUE)
  expect_match(a[1,2], "ref")
  
  expect_equivalent(attr(a, "rgroup"),
                    c("x", ""), info="Rgroup test")
  
  tmp <- getCrudeAndAdjustedModelData(fit)
  b <- printCrudeAndAdjustedModel(tmp, add_references=TRUE)
  
  expect_equivalent(a, b)
  # Getting the name wrong should not change the reference
  a <- printCrudeAndAdjustedModel(fit, add_references=TRUE, 
                                  add_references_pos=list(a=3))
  expect_match(a[1,2], "ref")
  
  # This should move the reference
  a <- printCrudeAndAdjustedModel(fit, add_references=TRUE, add_references_pos=list(x=2))
  expect_match(a[2,2], "ref")
  
  # Should end up at first position if referenced outside
  expect_warning(a <- printCrudeAndAdjustedModel(fit, add_references=TRUE, add_references_pos=list(x=5)))
  expect_match(a[1,2], "ref")
  
  # Bug with the same label occurring miultiple times
  complx_fit <- update(fit, .~x + boolean +
                         same_label +
                         same_labell +
                         same_labelll)
  a <- printCrudeAndAdjustedModel(complx_fit, 
                                  add_references=TRUE, desc_column=TRUE)
  expect_equal(nrow(a), 11)

  a2 <- printCrudeAndAdjustedModel(complx_fit)
  expect_lt(nrow(a2), nrow(a))
  
  a3 <- printCrudeAndAdjustedModel(complx_fit, desc_column = TRUE)
  expect_equivalent(nrow(a3), nrow(a), info = "When descriptive column is used then references should be added by default")
  
  a2 <- printCrudeAndAdjustedModel(complx_fit, desc_column = TRUE, order = c("same_label+","x", "boolean"))
  expect_equivalent(nrow(a2), nrow(a))
  expect_equivalent(tail(rownames(a2), 1), "boolean")
})

test_that("Test rbind", {
  set.seed(10)
  n <- 500
  ds <- data.frame(
    ftime = rexp(n),
    fstatus = sample(0:1, size = n, replace = TRUE),
    x = factor(sample(LETTERS[1:4], size = n, replace = TRUE)),
    same_label = factor(sample(c("Yes", "No"), size = n, replace = TRUE)),
    same_labell = factor(sample(c("Yes", "No"), size = n, replace = TRUE)),
    same_labelll = factor(sample(c("Yes", "No"), size = n, replace = TRUE)),
    boolean = sample(c(TRUE, FALSE), size = n, replace = TRUE),
    subsetting = factor(sample(c(TRUE, FALSE), size = n, replace = TRUE)))
  
  library(survival)
  fit1 <- coxph(Surv(ftime, fstatus == 1) ~ x + boolean, data=ds)
  fit2 <- coxph(Surv(ftime, fstatus == 1) ~ x + same_label, data=ds)
  
  a1 <- printCrudeAndAdjustedModel(fit1, add_references=TRUE)
  a2 <- printCrudeAndAdjustedModel(fit2, add_references=TRUE)
  a3 <- rbind(a1, a2)
  expect_equivalent(nrow(a3), nrow(a1) + nrow(a2))
  expect_true(is.null(attr(a3, "tspanner")))
  
  a3 <- rbind(a1 = a1, a2 = a2)
  expect_equivalent(nrow(a3), nrow(a1) + nrow(a2))
  expect_equivalent(attr(a3, "tspanner"), c("a1", "a2"))
  sink(file=ifelse(Sys.info()["sysname"] == "Windows",
                   "NUL",
                   "/dev/null"))
  expect_true(inherits(print(a3), "htmlTable"))
  sink()
})

test_that("Variable select",{
  set.seed(10)
  n <- 500
  ds <- data.frame(
    y = sample(0:1, size = n, replace = TRUE),
    x1 = factor(sample(LETTERS[1:4], size = n, replace = TRUE)),
    x2 = factor(sample(c("Yes", "No"), size = n, replace = TRUE)),
    x3 = factor(sample(c("Yes", "No"), size = n, replace = TRUE)),
    subsetting = factor(sample(c(TRUE, FALSE), size = n, replace = TRUE)))
  
  library(rms)
  dd <<- datadist(ds)
  options(datadist="dd")
  
  fit <- Glm(y ~ x1 + x2 + x3, data=ds, family=binomial)
  
  a <- printCrudeAndAdjustedModel(fit, order = c("x[12]"), add_references=TRUE)
  expect_equivalent(attr(a, "rgroup"), c("x1", "x2"))
  
  a <- printCrudeAndAdjustedModel(fit, order = c("x2", "x1"), add_references=TRUE)
  expect_equivalent(attr(a, "rgroup"), c("x2", "x1"))
  
  fit <- glm(y ~ x1 + x2 + x3, data=ds, family=binomial)
  
  a <- printCrudeAndAdjustedModel(fit, order = c("x[12]"), add_references=TRUE)
  expect_equivalent(attr(a, "rgroup"), c("x1", "x2"))
  
  a <- printCrudeAndAdjustedModel(fit, order = c("x2", "x1"), add_references=TRUE)
  expect_equivalent(attr(a, "rgroup"), c("x2", "x1"))
})

test_that("Check statistics",{
  set.seed(10)
  n <- 500
  ds <- data.frame(
    y = sample(0:1, size = n, replace = TRUE),
    x1 = factor(sample(LETTERS[1:4], size = n, replace = TRUE)),
    x2 = rnorm(n),
    subsetting = factor(sample(c(TRUE, FALSE), size = n, replace = TRUE)))
  ds$x1[sample(1:nrow(ds), size = 100)] <- NA
  ds$x2[sample(1:nrow(ds), size = 100)] <- NA
  
  fit <- glm(y~x1 + x2, data = ds)
  out <- printCrudeAndAdjustedModel(fit, desc_column = TRUE, 
                                    desc_args = caDescribeOpts(digits = 2))
  expect_true(out["A","Total"] == as.character(sum(ds$x1 == "A", na.rm=TRUE)))
  expect_match(out["x2","Total"], sprintf("%.2f", mean(ds$x2, na.rm=TRUE)))
  #TODO

  library(rms)
  dd <<- datadist(ds)
  options(datadist="dd")
  
  fit <- lrm(y~x1 + x2, data = ds)
  out <- printCrudeAndAdjustedModel(fit, desc_column = TRUE, 
                                    desc_args = caDescribeOpts(digits = 2))
  expect_true(out["A","Total"] == as.character(sum(ds$x1 == "A", na.rm=TRUE)))
  expect_match(out["x2","Total"], sprintf("%.2f", mean(ds$x2, na.rm=TRUE)))
})

test_that("Issue #5", {
  set.seed(1)
  data <- data.frame(outcome=rnorm(100),
                     sex=sample(c("Male","Female"),100,TRUE),
                     country=sample(c("USA","UK","AUS"),100,TRUE))
  
  fit <- lm(outcome ~ sex + country, data=data)
  out <- printCrudeAndAdjustedModel(fit,desc_column=TRUE)
  expect_equal(as.integer(out["Female","Total"]), 
               sum(data$sex == "Female"))
  expect_equal(as.integer(out["Male","Total"]), 
               sum(data$sex == "Male"))
  data$sex[1] <- NA
  out <- printCrudeAndAdjustedModel(fit,desc_column=TRUE)
  expect_equal(as.integer(out["Female","Total"]), 
               sum(data$sex == "Female", na.rm=TRUE))
  expect_equal(as.integer(out["Male","Total"]), 
               sum(data$sex == "Male", na.rm=TRUE))
})

test_that("Subsetting and bindings",{
  set.seed(10)
  n <- 500
  ds <- data.frame(
    y = sample(0:1, size = n, replace = TRUE),
    x1 = factor(sample(LETTERS[1:4], size = n, replace = TRUE)),
    x2 = rnorm(n),
    subsetting = factor(sample(c(TRUE, FALSE), size = n, replace = TRUE)))
  ds$x1[sample(1:nrow(ds), size = 100)] <- NA
  ds$x2[sample(1:nrow(ds), size = 100)] <- NA
  
  fit <- glm(y~x1 + x2, data = ds)
  out <- printCrudeAndAdjustedModel(fit, desc_column = TRUE, 
                                    desc_args = caDescribeOpts(digits = 2))
  expect_equivalent(dim(out[1,1:2]),
                    c(1,2))

  expect_equivalent(dim(out[,4:5]),
                    c(nrow(out),2))
  
  expect_equivalent(dim(out[3:4,]),
                    c(2, ncol(out)))

  ret <- cbind(out[,2:5], out[,1])  
  expect_equal(ncol(ret), 5)
})


test_that("Errors for printCrudeAndJust",{
  set.seed(10)
  n <- 500
  ds <- data.frame(
    y = sample(0:1, size = n, replace = TRUE),
    x1 = factor(sample(LETTERS[1:4], size = n, replace = TRUE)),
    x2 = rnorm(n),
    subsetting = factor(sample(c(TRUE, FALSE), size = n, replace = TRUE)))
  ds$x1[sample(1:nrow(ds), size = 100)] <- NA
  ds$x2[sample(1:nrow(ds), size = 100)] <- NA
  
  fit <- glm(y~x1 + x2, data = ds)
  expect_error(printCrudeAndAdjustedModel(NULL, desc_column = TRUE, 
                                          desc_args = caDescribeOpts(digits = 2)))
  expect_error(printCrudeAndAdjustedModel(NULL, desc_column = TRUE, 
                                          desc_args = "wrong argument"))
  
})
  