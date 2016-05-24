############################################################################
#     MLwiN MCMC Manual
#
# 22  Using the Structured MVN framework for models . . . . . . . . . . .341
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

# 22.1 MCMC theory for Structured MVN models . . . . . . . . . . . . . . 341

# 22.2 Using the SMVN framework in practice . . . . . . . . . . . . . . .344

library(R2MLwiN)
# MLwiN folder
mlwin <- getOption("MLwiN_path")
while (!file.access(mlwin, mode = 1) == 0) {
  cat("Please specify the root MLwiN folder or the full path to the MLwiN executable:\n")
  mlwin <- scan(what = character(0), sep = "\n")
  mlwin <- gsub("\\", "/", mlwin, fixed = TRUE)
}
options(MLwiN_path = mlwin)

# User's input if necessary

## Read tutorial data
data(tutorial, package = "R2MLwiN")

## Define the model

(mymodel <- runMLwiN(normexam ~ 1 + (1 | school) + (1 | student), estoptions = list(EstM = 1), data = tutorial))


## Structured MVN

(mymodel <- runMLwiN(normexam ~ 1 + (1 | school) + (1 | student), estoptions = list(EstM = 1, mcmcOptions = list(smvn = 1)), 
  data = tutorial))

# 22.3 Model Comparison and structured MVN models . . . . . . . . . . . .349

## Define the model

## Gibbs

(mymodel <- runMLwiN(normexam ~ 1 + standlrt + (1 | school) + (1 | student), estoptions = list(EstM = 1), data = tutorial))

## SMCMC

(mymodel <- runMLwiN(normexam ~ 1 + standlrt + (1 | school) + (1 | student), estoptions = list(EstM = 1, mcmcOptions = list(smcm = 1)), 
  data = tutorial))

## Structured MVN

(mymodel <- runMLwiN(normexam ~ 1 + standlrt + (1 | school) + (1 | student), estoptions = list(EstM = 1, mcmcOptions = list(smvn = 1)), 
  data = tutorial))

# 22.4 Assessing the need for the level 2 variance . . . . . . . . . . . 350

sixway(mymodel@chains[, "RP2_var_Intercept", drop = FALSE], "sigma2u0")

set.seed(1)
tutorial$temp <- double2singlePrecision(rnorm(4059))

## Define the model

## IGLS

(mymodel <- runMLwiN(temp ~ 1 + standlrt + (1 | school) + (1 | student), data = tutorial))

## Gibbs

(mymodel <- runMLwiN(temp ~ 1 + standlrt + (1 | school) + (1 | student), estoptions = list(EstM = 1), data = tutorial))

## Structured MVN

(mymodel <- runMLwiN(temp ~ 1 + standlrt + (1 | school) + (1 | student), estoptions = list(EstM = 1, mcmcOptions = list(smvn = 1)), 
  data = tutorial))

summary(mymodel@chains[, "RP2_var_Intercept"])
sixway(mymodel@chains[, "RP2_var_Intercept", drop = FALSE], "sigma2u0")

# Chapter learning outcomes . . . . . . . . . . . . . . . . . . . . . . .355





############################################################################
