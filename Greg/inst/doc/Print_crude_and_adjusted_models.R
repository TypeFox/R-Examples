## ---- message=FALSE------------------------------------------------------
library(datasets)
data(mtcars)
mtcars$am <- factor(mtcars$am, labels = c("Automatic", "Manual"))
fit <- lm(mpg ~ cyl + disp + hp + am, data = mtcars)
library(Greg)
printCrudeAndAdjustedModel(fit)

## ------------------------------------------------------------------------
printCrudeAndAdjustedModel(fit, 
                           digits = 1, 
                           add_references = TRUE,
                           rowname.fn = function(n){
  if (n == "disp")
    return("Displacement (cu.in.)")
  if (n == "hp")
    return("Gross horsepower")
  if (n == "cyl")
    return("No. cylinders")
  if (n == "am")
    return("Transmission")
  return(n)
})

## ---- message=FALSE------------------------------------------------------
library(Hmisc)
label(mtcars$disp) <- "Displacement (cu.in)"
label(mtcars$cyl) <- "No. cylinders"
label(mtcars$hp) <- "Gross horsepower"
label(mtcars$am) <- "Transmission"

printCrudeAndAdjustedModel(fit, 
                           digits = 1, 
                           add_references = TRUE)

## ------------------------------------------------------------------------
fit_mpg <- lm(mpg ~ cyl + disp + hp + am, data = mtcars)
fit_weight <- lm(wt ~ cyl + disp + hp + am, data = mtcars)
p_mpg <- printCrudeAndAdjustedModel(fit_mpg, digits = 1, add_references = TRUE)
p_weight <- printCrudeAndAdjustedModel(fit_weight, digits = 1, add_references = TRUE)
rbind("Miles per gallon" = p_mpg, 
      "Weight (1000 lbs)" = p_weight)

cbind("Miles per gallon" = p_mpg, 
      "Weight (1000 lbs)" = p_weight)


## ------------------------------------------------------------------------
p_mpg[,1:2]

p_mpg[1:2,]

## ------------------------------------------------------------------------
library("survival")

set.seed(10)
n <- 500
ds <- data.frame(
  ftime = rexp(n),
  fstatus = sample(0:1, size = n, replace = TRUE),
  y = rnorm(n = n),
  x1 = factor(sample(LETTERS[1:4], size = n, replace = TRUE)),
  x2 = rnorm(n, mean = 3, 2),
  x3 = rnorm(n, mean = 3, 2),
  x4 = factor(sample(letters[1:3], size = n, replace = TRUE)))

library(survival)
library(splines)
fit <- coxph(Surv(ds$ftime, ds$fstatus == 1) ~ x1 + ns(x2, 4) + x3 + strata(x4), data=ds)

printCrudeAndAdjustedModel(fit, add_references = TRUE)

# Note that the crude is with the strata
a <- getCrudeAndAdjustedModelData(fit)
a["x3", "Crude"] == 
  exp(coef(coxph(Surv(ds$ftime, ds$fstatus == 1) ~ x3 + 
                   strata(x4), data=ds)))

