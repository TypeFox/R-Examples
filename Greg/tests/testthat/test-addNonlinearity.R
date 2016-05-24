library(testthat)
library(splines)

context("Check addNonlinearity")
test_that("Check regular glm", {
  n <- 100
  set.seed(123)
  nl_ds <- data.frame(x=sample(seq(from = 0,
                                   to = pi*3,
                                   length.out = n),
                               size=n,
                               replace=TRUE),
                      sex = factor(sample(c("Male", "Female"), 
                                          size = n, 
                                          replace = TRUE)))
  nl_ds$y <- 
    sin(nl_ds$x)*2 + 
    (nl_ds$sex == "Male")*1 +
    rnorm(n = n, mean = 0, sd = .5)
  
  l_ds <- nl_ds
  l_ds$y <- 
    nl_ds$x*2 + 
    (nl_ds$sex == "Male")*1 +
    rnorm(n = n, mean = 0, sd = .5)
  
  vals <- sapply(2:7, function(x) AIC(glm(sprintf("y ~ ns(x, %d) + sex",
                                                  x), data=nl_ds)))
  expect_equivalent(AIC(addNonlinearity(glm(y ~ x + sex, data = nl_ds), 
                                        min_fn = AIC,
                                        flex_param = 2:7,
                                        variable = "x", spline_fn = "ns",
                                        workers = FALSE)),
                    min(vals))
  
  
  expect_equivalent(addNonlinearity(glm(y ~ x + sex, data = l_ds), 
                                    min_fn = AIC,
                                    flex_param = 2:7,
                                    variable = "x", spline_fn = "ns",
                                    workers = FALSE),
                    glm(y ~ x + sex, data = l_ds))
  
})
    
test_that("Check regular lm", {
  n <- 100
  set.seed(123)
  nl_ds <- data.frame(x=sample(seq(from = 0,
                                   to = pi*3,
                                   length.out = n),
                               size=n,
                               replace=TRUE),
                      sex = factor(sample(c("Male", "Female"), 
                                          size = n, 
                                          replace = TRUE)))
  nl_ds$y <- 
    sin(nl_ds$x)*2 + 
    (nl_ds$sex == "Male")*1 +
    rnorm(n = n, mean = 0, sd = .5)
  
  l_ds <- nl_ds
  l_ds$y <- 
    nl_ds$x*2 + 
    (nl_ds$sex == "Male")*1 +
    rnorm(n = n, mean = 0, sd = .5)
  
  vals <- sapply(2:7, function(x) AIC(lm(sprintf("y ~ ns(x, %d) + sex",
                                                  x), data=nl_ds)))
  expect_equivalent(AIC(addNonlinearity(lm(y ~ x + sex, data = nl_ds), 
                                        min_fn = AIC,
                                        flex_param = 2:7,
                                        variable = "x", spline_fn = "ns",
                                        workers = FALSE)),
                    min(vals))

  expect_equivalent(addNonlinearity(lm(y ~ x + sex, data = l_ds), 
                                    min_fn = AIC,
                                    flex_param = 2:7,
                                    variable = "x", spline_fn = "ns",
                                    workers = FALSE),
                    lm(y ~ x + sex, data = l_ds))
})

test_that("That rms-functions work", {
  n <- 100
  set.seed(123)
  nl_ds <- data.frame(x=sample(seq(from = 0,
                                   to = pi*3,
                                   length.out = n),
                               size=n,
                               replace=TRUE),
                      sex = factor(sample(c("Male", "Female"), 
                                          size = n, 
                                          replace = TRUE)))
  nl_ds$y <- 
    sin(nl_ds$x)*2 + 
    (nl_ds$sex == "Male")*1 +
    rnorm(n = n, mean = 0, sd = .5)
  
  l_ds <- nl_ds
  l_ds$y <- 
    nl_ds$x*2 + 
    (nl_ds$sex == "Male")*1 +
    rnorm(n = n, mean = 0, sd = .5)
  

  library(rms)
  vals <- sapply(3:7, function(x) AIC(ols(as.formula(sprintf("y ~ rcs(x, %d) + sex", x)), 
                                          data=nl_ds)))
  expect_equivalent(AIC(addNonlinearity(ols(y ~ x + sex, data = nl_ds), 
                                        min_fn = AIC,
                                        flex_param = 3:7,
                                        variable = "x", 
                                        spline_fn = "rcs",
                                        workers = FALSE)),
                    min(vals))
  
  expect_error(AIC(addNonlinearity(ols(y ~ x + sex, data = nl_ds), 
                                   min_fn = AIC,
                                   flex_param = 3:7,
                                   variable = "x", 
                                   spline_fn = "ns",
                                   workers = FALSE)))  
  
  expect_equivalent(addNonlinearity(ols(y ~ x + sex, data = l_ds), 
                                    min_fn = AIC,
                                    flex_param = 3:7,
                                    variable = "x", 
                                    spline_fn = "rcs",
                                    workers = FALSE),
                    ols(y ~ x + sex, data = l_ds))
})
