############################################################################
#     MLwiN MCMC Manual
#
# 2   Single Level Normal Response Modelling . . . . . . . . . . . . . . .21
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

## Choose IGLS algoritm for estimation
(mymodel1 <- runMLwiN(normexam ~ 1 + standlrt + (1 | student), data = tutorial))

# 2.1 Running the Gibbs Sampler . . . . . . . . . . . . . . . . . . . . . 26

## Choose MCMC algoritm for estimation
(mymodel2 <- runMLwiN(normexam ~ 1 + standlrt + (1 | student), estoptions = list(EstM = 1), data = tutorial))

estimates <- mymodel2@chains
par(mfrow = c(2, 2))
plot(1:niter(estimates), estimates[, "deviance"], xlab = "iteration", ylab = expression(paste("Est. of deviance")), 
  type = "l")
plot(1:niter(estimates), estimates[, "FP_Intercept"], xlab = "iteration", ylab = expression(paste("Est. of ", beta[0])), 
  type = "l")
plot(1:niter(estimates), estimates[, "FP_standlrt"], xlab = "iteration", ylab = expression(paste("Est. of ", beta[1])), 
  type = "l")
plot(1:niter(estimates), estimates[, "RP1_var_Intercept"], xlab = "iteration", ylab = expression(paste("Est. of ", 
  sigma[e0]^2)), type = "l")

# 2.2 Deviance statistic and the DIC diagnostic . . . . . . . . . . . . . 28

# 2.3 Adding more predictors . . . . . . . . . . . . . . . . . . . . . . .29

## Choose IGLS algoritm for estimation
(mymodel3 <- runMLwiN(normexam ~ 1 + standlrt + sex + schgend + (1 | student), data = tutorial))

## Choose MCMC algoritm for estimation
(mymodel4 <- runMLwiN(normexam ~ 1 + standlrt + sex + schgend + (1 | student), estoptions = list(EstM = 1), data = tutorial))

# 2.4 Fitting school effects as fixed parameters . . . . . . . . . . . . .32

tutorial$school <- as.factor(tutorial$school)

## Choose MCMC algoritm for estimation (IGLS will be used to obtain starting values for MCMC)
(mymodel5 <- runMLwiN(normexam ~ 1 + standlrt + sex + school + (1 | student), estoptions = list(EstM = 1), data = tutorial))

# Chapter learning outcomes . . . . . . . . . . . . . . . . . . . . . . . 33





############################################################################
