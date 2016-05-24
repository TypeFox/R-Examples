context("plotHR private functions")

test_that("Check the prPhNewData function",{
  set.seed(10)
  n <- 11
  ds <- data.frame(
    ftime = rexp(n),
    fstatus = sample(0:1, size = n, replace = TRUE),
    x1 = factor(sample(c("a", "aa", "aaa"), size = n, replace = TRUE)),
    x2 = rnorm(n, mean = 3, 2),
    x3 = factor(sample(c("a", "aa", "aaa"), size = n, replace = TRUE)))

  library(rms)
  dd <<- datadist(ds)
  options(datadist = "dd")
  
  fit <- coxph(Surv(ftime, fstatus) ~ x1 + x2, data=ds)
  expect_error(prPhNewData(fit, "aa"))
  
  expect_equivalent(nrow(prPhNewData(fit, "x1")),
                    length(unique(ds$x1)),
                    info="The length should be equal to the different type within the data")

  expect_equivalent(nrow(prPhNewData(fit, "x2")),
                    length(unique(ds$x2)),
                    info="The length should be equal to the different type within the data")
  
  expect_equivalent(length(unique(prPhNewData(fit, "x2")$x1)),
                    1, 
                    info="All other variables should be only of one type")
  
  expect_equivalent(nrow(prPhNewData(fit, "x2", xlim=quantile(ds$x2, c(.05, .5)))),
                    length(unique(ds$x2[ds$x2 <= quantile(ds$x2, .5) &
                                          ds$x2 >= quantile(ds$x2, .05)])),
                    info="The values beyond the xlim should be ignored")
  
  fit <- cph(Surv(ftime, fstatus) ~ x1 + x2, data=ds)
  expect_error(prPhNewData(fit, "aa"))
  
  expect_equivalent(nrow(prPhNewData(fit, "x1")),
                    length(unique(ds$x1)),
                    info="The length should be equal to the different type within the data")
  
  expect_equivalent(nrow(prPhNewData(fit, "x2")),
                    length(unique(ds$x2)),
                    info="The length should be equal to the different type within the data")
  
  expect_equivalent(length(unique(prPhNewData(fit, "x2")$x1)),
                    1, 
                    info="All other variables should be only of one type")
  
  expect_equivalent(nrow(prPhNewData(fit, "x2", xlim=quantile(ds$x2, c(.05, .5)))),
                    length(unique(ds$x2[ds$x2 <= quantile(ds$x2, .5) &
                                          ds$x2 >= quantile(ds$x2, .05)])),
                    info="The values beyond the xlim should be ignored")
  
  S <- with(ds, Surv(ftime, fstatus))
  fit <- cph(S ~ x1 + x2, data=ds)
  expect_equivalent(ncol(prPhNewData(fit, "x2")),
                    2)
})



test_that("Check the prPhEstimate function",{
  # From the cph example
  library(rms)
  n <- 100
  set.seed(731)
  age <- 50 + 12*rnorm(n)
  label(age) <- "Age"
  sex <- factor(sample(c('Male','Female'), n, 
                       rep=TRUE, prob=c(.6, .4)))
  cens <- 15*runif(n)
  h <- .02*exp(.04*(age-50)+.8*(sex=='Female'))
  dt <- -log(runif(n))/h
  label(dt) <- 'Follow-up Time'
  e <- ifelse(dt <= cens,1,0)
  dt <- pmin(dt, cens)
  units(dt) <- "Year"
  dd <<- datadist(age, sex)
  options(datadist='dd')
  S <- Surv(dt,e)
  
  fit <- cph(S ~ age + sex, x=TRUE, y=TRUE)
  est <- prPhEstimate(fit, 
                      alpha=.05,
                      term.label = "age",
                      ylog=TRUE,
                      cntrst=TRUE)
  expect_equivalent(colnames(est),
                    c("xvalues", "estimate",
                      "lower", "upper"))
  
  est <- prPhEstimate(fit, 
                      alpha=.05,
                      term.label = "age",
                      ylog=TRUE,
                      cntrst=FALSE)
  expect_equivalent(colnames(est),
                    c("xvalues", "estimate",
                      "lower", "upper"))
  

  est <- prPhEstimate(fit, 
                      alpha=.05,
                      term.label = "age",
                      ylog=TRUE,
                      cntrst=TRUE)

  expect_true(coef(fit)["age"] *
                (with(est, estimate[which.max(xvalues)] - 
                        estimate[which.min(xvalues)])) > 0,
              info="Two negative estimates should produce a positive value as well as two positive")

  fit <- coxph(S ~ age + sex, x=TRUE, y=TRUE)
  expect_error(prPhEstimate(fit, 
                            term.label = "age",
                            alpha=.05,
                            ylog=TRUE,
                            cntrst=TRUE))

  est <- prPhEstimate(fit, 
                      term.label = "age",
                      alpha=.05,
                      ylog=TRUE,
                      cntrst=FALSE)
  expect_equivalent(colnames(est),
                    c("xvalues", "estimate",
                      "lower", "upper"))

  est <- prPhEstimate(fit, 
                      term.label = "age",
                      alpha=.05,
                      ylog=FALSE,
                      cntrst=FALSE)
  expect_true(all(est$estimate > 0))
  expect_true(all(est$lower > 0))
  expect_true(all(est$upper > 0))
  expect_true(all(est$xvalues <= max(age)))
  expect_true(all(est$xvalues >= min(age)))
  
  fit <- cph(S ~ rcs(age, 3) + sex, x=TRUE, y=TRUE)
  est <- prPhEstimate(fit, 
                      alpha=.05,
                      term.label = "age",
                      ylog=TRUE,
                      cntrst=TRUE)
  expect_equivalent(colnames(est),
                    c("xvalues", "estimate",
                      "lower", "upper"))
  
  fit <- coxph(S ~ pspline(age, 3) + sex)
  est <- prPhEstimate(fit, 
                      alpha=.05,
                      term.label = "age",
                      ylog=TRUE,
                      cntrst=FALSE)
  expect_equivalent(colnames(est),
                    c("xvalues", "estimate",
                      "lower", "upper"))
  
})
  
test_that("General plotting funcitonality of plotHR",{
  library(survival)
  library(rms)
  
  # Get data for example
  n <- 1000
  set.seed(731)
  
  age <- round(50 + 12*rnorm(n), 1)
  label(age) <- "Age"
  
  sex <- factor(sample(c('Male','Female'), n, 
                       rep=TRUE, prob=c(.6, .4)))
  cens <- 15*runif(n)
  
  smoking <- factor(sample(c('Yes','No'), n, 
                           rep=TRUE, prob=c(.2, .75)))
  
  h <- .02*exp(.02*(age-50)+.1*((age-50)/10)^3+.8*(sex=='Female')+2*(smoking=='Yes'))
  dt <- -log(runif(n))/h
  label(dt) <- 'Follow-up Time'
  
  e <- ifelse(dt <= cens,1,0)
  dt <- pmin(dt, cens)
  units(dt) <- "Year"
  
  # Add missing data to smoking
  smoking[sample(1:n, round(n*0.05))] <- NA
  
  # Create a data frame since plotHR will otherwise
  # have a hard time getting the names of the variables
  ds <- data.frame(
    dt = dt,
    e = e,
    age=age, 
    smoking=smoking,
    sex=sex)
  
  library(splines)
  Srv <- Surv(dt,e)
  fit.coxph <- coxph(Srv ~ bs(age, 3) + sex + smoking, data=ds)
  
  org_par <- par(xaxs="i")
  plotHR(fit.coxph, term="age", plot.bty="o", xlim=c(30, 70), xlab="Age")
  
  dd <- datadist(ds)
  options(datadist='dd')
  fit.cph <- cph(Srv ~ rcs(age,4) + sex + smoking, data=ds, x=TRUE, y=TRUE)
  
  plotHR(fit.cph, term=1, plot.bty="l", xlim=c(30, 70), xlab="Age")
  
  plotHR(fit.cph, term="age", plot.bty="l", xlim=c(30, 70), ylog=FALSE, rug="ticks", xlab="Age")
  
  unadjusted_fit <- cph(Srv ~ rcs(age,4), data=ds, x=TRUE, y=TRUE)
  plotHR(list(fit.cph, unadjusted_fit), term="age", xlab="Age",
         polygon_ci=c(TRUE, FALSE), 
         col.term = c("#08519C", "#77777799"),
         col.se = c("#DEEBF7BB", grey(0.6)),
         lty.term = c(1, 2),
         plot.bty="l", xlim=c(30, 70))
})