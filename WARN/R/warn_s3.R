# =====R code=========================
# Weaning age reconstruction using d15N.
# ==============================
# 2013-01-14: Tsutaya T: combined all files,
#  integrated names of variable.
# https://sites.google.com/site/leihcrev/r/writing-your-own-functions
# 2014-11-03: Tsutaya T: added stop() message for 'age' > 10.
# ==============================
# OBJECTIVE ----------
# This program performs Apporoximate Bayesian Computation with SMC
#  adopting Partial Rejection Control (PRC)
#  in order to calculates posterior distributuions
#  of the weaning ages and weaning parameters.
# ==============================
# FUNCTION DEFINITIONS ----------
# ==============================
# Define the class "warn".
# ==============================
warn <- function(age, d15N, female.mean, female.sd,
  prior, num.particle, form, tolerances){
# S3 method.
  UseMethod("warn")
}

warn.default <- function(age, d15N, female.mean, female.sd = NA,
  prior = c(0.5, 3, 3, 3, 1.9, 0.9, female.mean, 3),
  num.particle = 10000, form = "parabolic",
  tolerances = c(2.0, 1.0, 0.5, 0.25, 0.125, 0.0625, 0)){
# Default method for S3 "warn" class (= WARN()).
#
# arg:
#  all arguments will be passed to WARN().
#
# return:
#  Class "warn" object.
  if(max(age) > 10){
    stop(message=">10 years of 'age'")
  }

  warn.result <- WARN(age = age,
    d15N = d15N,
    female.mean = female.mean,
    female.sd = female.sd,
    prior = prior,
    num.particle = num.particle,
    form = form,
    tolerances = tolerances)
  warn.result$call <- match.call()
  class(warn.result) <- "warn"

  return(warn.result)
}

print.warn <- function(x, ...){
# Print method for S3 "warn" class.
#
# arg:
#  x: An object which was defined as "warn" class.
#
# returns:
#  Pinrt of MDEs for the object.
  cat("Call:\n")
  print(x$call) # Which method did match.
  cat("\nMDEs:\n")
  print(x$mde)
  cat("\nJoint probability of the waening ages:\n")
  print(x$prob.2d.age)

  cat("\n")
}

plot.warn <- function(x, 
  hline.female = TRUE, hline.adult = FALSE, 
  adult.mean = NA, adult.sd = 0,
  is.legend = TRUE, is.female = TRUE, ...){
# Plot method for S3 "warn" class.
#
# arg:
#  x: An object which was defined as "warn" class.
#  ...: Passed to plot() in DrawMDE.
#   i.e. hline, adult.mean, adult.sd, is.legend, is.female, ...
#
# function:
#  DrawMDE
#
# returns:
#  plot of summaries for the object.
  DrawMDE(par.mde = x$mde[ , 1],
    d15N = x$d15N,
    age = x$age,
    female.mean = x$female.mean,
    female.sd = x$female.sd,
    form = x$form,
    hline.female = hline.female,
    hline.adult = hline.adult, 
    adult.mean = adult.mean,
    adult.sd = adult.sd,
    is.legend = is.legend,
    is.female = is.female, ...)
}

summary.warn <- function(object, ...){
# Summary method for S3 "warn" class.
#
# arg:
#  object: An object which was defined as "warn" class.
#  \dots: Additional arguments affecting the summary produced.
#
# returns:
#  Print of class "warnSummary" object.
  object$individual <- length(object$age)
  class(object) <- "warnSummary"
  return(object)
}
# ==============================
# Define the class "warnSummary".
# ==============================
warnSummary <- function(object, ...){
  UseMethod("warnSummary")
}

warnSummary.default <- function(object, ...){
# Default method for S3 "warnSummary" class (= summary(warn)).
#
# arg:
#  object: an object of class "warnSummary", will be passed to print.warnSummary().
#  \dots: Additional arguments affecting the summary produced.
#
# return:
#  Summary of class "warnSummary" object.
  result <- summary(object = object)
  return(result)
}

print.warnSummary <- function(x, ...){
# Print method for S3 "warnSummary" class.
#
# arg:
#  x: An object which was defined as "warnSummary" class.
#
# returns:
#  Pinrt of summaries for the object.
  cat("Call:\n")
  print(x$call) # Which method did match.
  cat("\nMDEs and marginal probabilities:\n")
  print(x$mde)
  cat("\nJoint probability of the waening ages:\n")
  print(x$prob.2d.age)
  cat("\nMean squared distance under MDE parameters:\n")
  print(x$dist.mde)
  cat("\nNumber of non-adult individuals:\n")
  print(x$individual)
  cat("\nNumber of particles:\n")
  print(x$particle)

  cat("\n")
}

# ==============================
# Define the class "warnProb".
# ==============================
warnProb <- function(object, weaning.par, range.x, range.y){
  UseMethod("warnProb")
}

warnProb.default <- function(object, weaning.par = "age",
  range.x, range.y = NA){
# Default method for S3 "warnProb" class (= CalcProb).
#
# arg:
#  object: An object which was defined as "warnProb" class.
#  weaning.par: Weaning parameters interested, c("age", "enrich", "wnfood").
#  range: A range with which we want to calculate probability.
#   Fractional point lower than e-002 is omitted.
#
# return:
#  Class "warnProb" object.
  # Calculate probabilities.
  if(weaning.par == "age"){
    warnprob.result <- CalcProb2D(
      kde = object$kde.age,
      range.x = range.x,
      range.y = range.y)
  }else if(weaning.par == "enrich"){
    warnprob.result <- CalcProb1D(
      kde = object$kde.enrich,
      range.x = range.x)
    if(!is.na(range.y[1])){
      warning("Argument for range.y is ignored.")
    }
  }else if(weaning.par == "wnfood"){
    warnprob.result <- CalcProb1D(
      kde = object$kde.wnfood,
      range.x = range.x)
    if(!is.na(range.y[1])){
      warning("Argument for range.y is ignored.")
    }
  }else{
    stop("Choose correct weaning paramters from c(\"age\", \"enrich\", \"wnfood\").\n")
  }

  # Another slots.
  warnprob.result <- c(object, warnprob.result)
  warnprob.result$weaning.par <- weaning.par
  warnprob.result$call <- match.call()
  class(warnprob.result) <- "warnProb"

  return(warnprob.result)
}

print.warnProb <- function(x, ...){
# Print method for S3 "warnProb" class.
#
# arg:
#  x: An object which was defined as "warnProb" class.
#
# returns:
#  Print of probability for the selected ranges.
  cat("Call:\n")
  print(x$call) # Which method did match.
  cat("\nProbability:\n")
  print(x$probability)

  cat("\n")
}

plot.warnProb <- function(x,
  is.legend = TRUE, is.contour = TRUE,
  is.image = FALSE, is.prior = FALSE, ...){
# Plot method for S3 "warnProb" class.
#
# arg:
#  x: An object which was defined as "warnProb" class.
#  ...: Passed to DrawProb (= plot).
#   i.e. age: is.legend, is.image, ...
#   i.e. enrich, wnfood: is.legend, is.prior, ...
#
# function:
#  DrawProb1D
#  DrawProb2D
#
# returns:
#  plot of summaries for the probability estimation.
#
# note:
#  Dispatch for weaning paramters is working.
  if(x$weaning.par == "age"){
    DrawProb2D(kde = x$kde,
      range.x1 = x$range[1:2],
      range.y1 = x$range[3:4],
      mde = x$mde[1:2, 1],
      is.legend = is.legend,
      is.contour = is.contour,
      is.image = is.image, ...)
  }else if(x$weaning.par == "enrich"){
    DrawProb1D(kde = x$kde,
      range.x1 = x$range,
      mde = x$mde[3, 1],
      is.legend = is.legend,
      is.prior = is.prior,
      hyper.par = x$prior[5:6], ...)
  }else if(x$weaning.par == "wnfood"){
    DrawProb1D(kde = x$kde,
      range.x1 = x$range,
      mde = x$mde[4, 1],
      is.legend = is.legend,
      is.prior = is.prior,
      hyper.par = x$prior[7:8], ...)
  }else{
    stop("Weaning paramter is incorrect.\n")
  }
}

summary.warnProb <- function(object, ...){
# Summary method for S3 "warnProb" class.
#
# arg:
#  object: An object which was defined as "warnProb" class.
#  \dots: Additional arguments affecting the summary produced.
#
# returns:
#  Print of class "warnProbSummary" object.
  class(object) <- "warnProbSummary"
  return(object)
}

# ==============================
# Define the class "warnSummary".
# ==============================
warnProbSummary <- function(object, ...){
  UseMethod("warnProbSummary")
}

warnProbSummary.default <- function(object, ...){
# Default method for S3 "warnProbSummary" class (= summary(warnProb)).
#  \dots: Additional arguments affecting the summary produced.
#
# arg:
#  object: an object of class "warnProbSummary",
#   will be passed to print.warnProbSummary().
#
# return:
#  Summary of class "warnSummary" object.
  result <- summary(object = object, ...)
  return(result)
}

print.warnProbSummary <- function(x, ...){
# Print method for S3 "warnSummary" class.
#
# arg:
#  x: An object which was defined as "warnProbSummary" class.
#
# returns:
#  Pinrt of summaries for the object.
  cat("Call:\n")
  print(x$call) # Which method did match.
  cat("\nWeaning parameter:\n")
  print(x$weaning.par)
  cat("\nRange:\n")
  print(x$range)
  cat("\nProbability:\n")
  print(x$probability)

  cat("\n")
}

# ==============================
# Define the class "warnOptim".
# ==============================
warnOptim <- function(age, ...){
# S3 method.
  UseMethod("warnOptim")
}

warnOptim.default <- function(age, d15N, female.mean,
  par.initial = c(0.5, 3, 1.9, female.mean),
  form = "parabolic", ...){
# Default method for S3 "warnOptim" class (= WrapperOptim).
#
# arg:
#  all arguments will be passed to WrapperOptim().
#
# return:
#  Class "warnOptim" object.
  warnoptim.result <- WrapperOptim(age = age,
    d15N = d15N,
    female.mean = female.mean,
    par.initial = par.initial,
    form = form, ...)
  warnoptim.result$call <- match.call()
  class(warnoptim.result) <- "warnOptim"

  return(warnoptim.result)
}

print.warnOptim <- function(x, ...){
# Print method for S3 "warnOptim" class.
#
# arg:
#  x: An object which was defined as "warnOptim" class.
#
# function:
#  ConditionAge
#
# returns:
#  Pinrt of values of the optim function.
	# Conditioning the weaning ages.
	weaning.age <- ConditionAge(
    t1 = x$par[1], t2 = x$par[2])
	x$par[1] <- weaning.age[1]
	x$par[2] <- weaning.age[2]

  par.named <- x$par
  names(par.named) <- c("t1", "t2", "enrich", "wnfood")

  cat("Call:\n")
  print(x$call) # Which method did match.
  cat("\nOptimized parameters:\n")
  print(par.named)

  cat("\n")
}

# ==============================
# Define the class "warnCI".
# ==============================
warnCI <- function(object, threshold){
  UseMethod("warnCI")
}

warnCI.default <- function(object, threshold = 0.95){
# Default method for S3 "warnCI" class (= CalcCI).
#
# arg:
#  object: An object which was defined as "warnCI" class.
#  threshold: A scalar or vector indicating threshold of CI
#   for age, enrich, wnfood. If scalar, value is repeated.
#   From 0 to 1.
#
# return:
#  Class "warnCI" object.
  # Threshold.
  if(length(threshold) == 1){
    threshold <- rep(threshold, 3)
  }

  # Calculate CIs.
  warnci.age <- CalcCI2D(
    kde = object$kde.age,
    mde.x = object$mde[1, 1],
    mde.y = object$mde[2, 1],
    threshold = threshold[1])
  warnci.enrich <- CalcCI1D(
    kde = object$kde.enrich,
    mde.x = object$mde[3, 1],
    threshold = threshold[2])
  warnci.wnfood <- CalcCI1D(
    kde = object$kde.wnfood,
    mde.x = object$mde[4, 1],
    threshold = threshold[3])
  names(threshold) <- c("age", "enrich", "wnfood")

  warnci.result <- list(
    ci.age = warnci.age,
    ci.enrich = warnci.enrich,
    ci.wnfood = warnci.wnfood,
    ci.threshold = threshold)
  
  # Another slots.
  result <- c(object, warnci.result)
  result$call <- match.call()
  class(result) <- "warnCI"

  return(result)
}

print.warnCI <- function(x, ...){
# Print method for S3 "warnCI" class.
#
# arg:
#  x: An object which was defined as "warnCI" class.
#
# returns:
#  Print of calculated ranges and its probability.
  cat("Call:\n")
  print(x$call) # Which method did match.
  cat("\nCI and its probability:\n")
  cat("\nWeaning ages:\n")
  print(x$ci.age)
  cat("\nEnrichment of d15N from mother:\n")
  print(x$ci.enrich)
  cat("\nd15N of collagen derived entirelly from weaning foods:\n")
  print(x$ci.wnfood)

  cat("\n")
}

plot.warnCI <- function(x, weaning.par = "age", ...){
# Plot method for S3 "warnCI" class.
#
# arg:
#  x: An object which was defined as "warnCI" class.
#  weaning.par: Weaning parameters interested, c("age", "enrich", "wnfood").
#  \dots: Passed to DrawCI (= plot).
#   i.e. age: is.legend, is.image, ...
#   i.e. enrich, wnfood: is.legend, is.prior, ...
#
# function:
#  DrawCI1D
#  DrawCI2D
#
# returns:
#  plot of summaries for the probability estimation.
#
# note:
#  Dispatch for weaning paramters is working.
  if(weaning.par == "age"){
    DrawProb2D(kde = x$kde.age,
      range.x1 = x$ci.age$range[1:2],
      range.y1 = x$ci.age$range[3:4],
      mde = x$mde[1:2, 1], ...)
  }else if(weaning.par == "enrich"){
    DrawProb1D(kde = x$kde.enrich,
      range.x1 = x$ci.enrich$range,
      mde = x$mde[3, 1], ...)
  }else if(weaning.par == "wnfood"){
    DrawProb1D(kde = x$kde.wnfood,
      range.x1 = x$ci.wnfood$range,
      mde = x$mde[4, 1], ...)
  }else{
    stop("Choose correct weaning paramters from c(\"age\", \"enrich\", \"wnfood\").\n")
  }
}

