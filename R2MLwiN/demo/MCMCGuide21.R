############################################################################
#     MLwiN MCMC Manual
#
# 21  Using Structured MCMC . . . . . . . . . . . . . . . . . . . . . . .327
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

# 21.1 SMCMC Theory . . . . . . . . . . . . . . . . . . . . . . . . . . .327

# 21.2 Fitting the model using MLwiN . . . . . . . . . . . . . . . . . . 330

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

summary(mymodel@chains[, "FP_Intercept"])
sixway(mymodel@chains[, "FP_Intercept", drop = FALSE], "beta_0")

## Structured MCMC

(mymodel <- runMLwiN(normexam ~ 1 + (1 | school) + (1 | student), estoptions = list(EstM = 1, mcmcOptions = list(smcm = 1)), 
  data = tutorial))

summary(mymodel@chains[, "FP_Intercept"])
sixway(mymodel@chains[, "FP_Intercept", drop = FALSE], "beta_0")

# 21.3 A random intercepts model . . . . . . . . . . . . . . . . . . . . 334

(mymodel <- runMLwiN(normexam ~ 1 + standlrt + (1 | school) + (1 | student), estoptions = list(EstM = 1, mcmcOptions = list(smcm = 1)), 
  data = tutorial))

trajectories(mymodel, Range = c(1, 500))

# 21.4 Examining the residual chains . . . . . . . . . . . . . . . . . . 335

(mymodel <- runMLwiN(normexam ~ 1 + standlrt + (1 | school) + (1 | student), estoptions = list(EstM = 1, resi.store = TRUE, 
  resi.store.levs = 2, mcmcMeth = list(iterations = 5001), mcmcOptions = list(smcm = 1)), data = tutorial))

## Each row represents each iteration
sixway(mymodel@resi.chains$resi_lev2[, 1, drop = FALSE], name = "school1")

# 21.5 Random slopes model theory . . . . . . . . . . . . . . . . . . . .336

# 21.6 Random Slopes model practice . . . . . . . . . . . . . . . . . . .338

(mymodel <- runMLwiN(normexam ~ 1 + standlrt + (1 + standlrt | school) + (1 | student), estoptions = list(EstM = 1, 
  mcmcOptions = list(smcm = 1)), data = tutorial))

sixway(mymodel@chains[, "FP_Intercept", drop = FALSE], "beta_0")
sixway(mymodel@chains[, "FP_standlrt", drop = FALSE], "beta_1")

# Chapter learning outcomes . . . . . . . . . . . . . . . . . . . . . . .340





############################################################################
