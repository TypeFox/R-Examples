library("testthat")
context("forestplotRegrObj")

# simulated data to test 
set.seed(1000)
n <- 1000
cov <- data.frame(
  ftime = rexp(n),
  fstatus = sample(0:2,n,replace=TRUE),
  x1 = runif(n),
  x2 = runif(n),
  x3 = runif(n))

library(rms)
dd <<- datadist(cov)
options(datadist="dd")

test_that("Basic test for coverage for forestplotRegrObj", {
  # TODO: Add more specific tests
  fit1 <- cph(Surv(ftime, fstatus == 1) ~ x1 + x2 + x3, data=cov)
  fit2 <- cph(Surv(ftime, fstatus == 2) ~ x1 + x2 + x3, data=cov)

  forestplotRegrObj (regr.obj = fit1, new_page=TRUE)
  
  library(forestplot)
  forestplotRegrObj (regr.obj = list(fit1, fit2),
                     legend = c("Status = 1", "Status = 2"), 
                     legend_args = fpLegend(title="Type of regression"),
                     new_page=TRUE)
  
  modifyNameFunction <- function(x){
    if (x == "x1")
      return ("Covariate A")
    
    if (x == "x2")
      return (expression(paste("My ", beta[2])))
    
    return (x)
  }
  
  forestplotRegrObj (regr.obj = list(fit1, fit2), 
                     col=fpColors(box=c("darkblue", "darkred")),
                     variablesOfInterest.regexp = "(x2|x3)",
                     legend = c("First model", "Second model"),
                     legend_args = fpLegend(title = "Models"),
                     rowname.fn = modifyNameFunction, new_page=TRUE)
  
  forestplotRegrObj (regr.obj = list(fit1, fit2), 
                     col=fpColors(box=c("darkblue", "darkred")),
                     variablesOfInterest.regexp = "(x2|x3)",
                     order.regexps = c("x3", "x2"),
                     legend = c("First model", "Second model"),
                     legend_args = fpLegend(title = "Models"),
                     rowname.fn = modifyNameFunction, new_page=TRUE)
})
