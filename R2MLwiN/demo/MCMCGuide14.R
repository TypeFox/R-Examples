############################################################################
#     MLwiN MCMC Manual
#
# 14  Adjusting for Measurement Errors in Predictor Variables . . . . . .199
#
#     Browne, W.J. (2009) MCMC Estimation in MLwiN, v2.13. Centre for
#     Multilevel Modelling, University of Bristol.
############################################################################
#     R script to replicate all analyses using R2MLwiN
#
#     Zhang, Z., Charlton, C., Parker, R, Leckie, G., and Browne, W.J.
#     Centre for Multilevel Modelling, 2012
#     http://www.bristol.ac.uk/cmm/software/R2MLwiN/
############################################################################

library(R2MLwiN)
# MLwiN folder
mlwin <- getOption("MLwiN_path")
while (!file.access(mlwin, mode = 1) == 0) {
  cat("Please specify the root MLwiN folder or the full path to the MLwiN executable:\n")
  mlwin <- scan(what = character(0), sep = "\n")
  mlwin <- gsub("\\", "/", mlwin, fixed = TRUE)
}
options(MLwiN_path = mlwin)

## Read tutorial data
data(tutorial, package = "R2MLwiN")

# 14.1 Effects of measurement error on predictors . . . . . . . . . . . .200
set.seed(1)
error <- double2singlePrecision(rnorm(length(tutorial$standlrt), 0, sqrt(0.2)))
obslrt <- double2singlePrecision(tutorial$standlrt + error)
tutorial <- cbind(tutorial, error, obslrt)

(mymodel <- runMLwiN(normexam ~ 1 + standlrt + (1 | student), data = tutorial))

(mymodel <- runMLwiN(normexam ~ 1 + error + (1 | student), data = tutorial))

(mymodel <- runMLwiN(normexam ~ 1 + obslrt + (1 | student), data = tutorial))

(mymodel <- runMLwiN(normexam ~ 1 + obslrt + (1 | student), estoptions = list(EstM = 1, merr = c(N = 1, "obslrt",
  0.2)), data = tutorial))

# 14.2 Measurement error modelling in multilevel models . . . . . . . . .205

(mymodel <- runMLwiN(normexam ~ 1 + standlrt + (1 + standlrt | school) + (1 | student), estoptions = list(EstM = 1),
  data = tutorial))

(mymodel <- runMLwiN(normexam ~ 1 + obslrt + (1 + obslrt | school) + (1 | student), estoptions = list(EstM = 1), data = tutorial))

(mymodel <- runMLwiN(normexam ~ 1 + obslrt + (1 + obslrt | school) + (1 | student), estoptions = list(EstM = 1, merr = c(N = 1,
  "obslrt", 0.2)), data = tutorial))

# 14.3 Measurement errors in binomial models . . . . . . . . . . . . . . 208


## Read bang1 data
data(bang1, package = "R2MLwiN")

set.seed(1)
bang1$obsage <- double2singlePrecision(bang1$age + rnorm(length(bang1$age), 0, 5))

(mymodel <- runMLwiN(logit(use) ~ 1 + age, D = "Binomial", estoptions = list(EstM = 1), data = bang1))

(mymodel <- runMLwiN(logit(use) ~ 1 + obsage, D = "Binomial", estoptions = list(EstM = 1), data = bang1))

## Adjust for the measurement errors
(mymodel <- runMLwiN(logit(use) ~ 1 + obsage, D = "Binomial", estoptions = list(EstM = 1, merr = c(N = 1,
  "obsage", 25)), data = bang1))

# 14.4 Measurement errors in more than one variable and
#      misclassifications . . . . . . . . . . . . . . . . . . . . . . . .211

# Chapter learning outcomes . . . . . . . . . . . . . . . . . . . . . . .128





############################################################################
