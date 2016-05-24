library("testthat")
context("simpleRmsAnova")

# ** Borrowed code from the lrm example **
library(rms)
n <- 1000    # define sample size
set.seed(17) # so can reproduce the results
age            <- rnorm(n, 50, 10)
blood.pressure <- rnorm(n, 120, 15)
cholesterol    <- rnorm(n, 200, 25)
sex            <- factor(sample(c('female','male'), n,TRUE))
label(age)            <- 'Age'      # label is in Hmisc
label(cholesterol)    <- 'Total Cholesterol'
label(blood.pressure) <- 'Systolic Blood Pressure'
label(sex)            <- 'Sex'
units(cholesterol)    <- 'mg/dl'   # uses units.default in Hmisc
units(blood.pressure) <- 'mmHg'

#To use prop. odds model, avoid using a huge number of intercepts by
#grouping cholesterol into 40-tiles

# Specify population model for log odds that Y=1
L <- .4*(sex=='male') + .045*(age-50) +
  (log(cholesterol - 10)-5.2)*(-2*(sex=='female') + 2*(sex=='male'))
# Simulate binary y to have Prob(y=1) = 1/[1+exp(-L)]
y <- ifelse(runif(n) < plogis(L), 1, 0)
cholesterol[1:3] <- NA   # 3 missings, at random

ddist <- datadist(age, blood.pressure, cholesterol, sex)
options(datadist='ddist')

test_that("Basic test for coverage for simpleRmsAnova", {
  # TODO: Add more specific tests
  
  fit_lrm <- lrm(y ~ blood.pressure + sex * (age + rcs(cholesterol,4)),
                 x=TRUE, y=TRUE)
  
  a_out <- anova(fit_lrm, 
                 dec.F = 1,
                 ss    = FALSE)
  
  a <- simpleRmsAnova(a_out, 
                      subregexps = rbind(c("blood.pressure", "Blood pressure"),
                                         c("age", "Age"),
                                         c("cholesterol", "Cholesterol"),
                                         c("sex", "Sex")),
                      caption="Anova output for a logistic regression model")
  expect_equivalent(dim(a), dim(a_out))
  expect_equivalent(sum(grepl("Blood pressure", dimnames(a)[[1]], fixed = TRUE)), 1)

  sink(file=ifelse(Sys.info()["sysname"] == "Windows",
                   "NUL",
                   "/dev/null"))
  p_out <- print(a)
  sink()
  expect_true(inherits(p_out, "htmlTable"))
  
  a <- simpleRmsAnova(a_out, 
                      caption="Anova output for a logistic regression model")
  expect_equivalent(sum(grepl("Blood pressure", dimnames(a)[[1]], fixed = TRUE)), 0)
})
