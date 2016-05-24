############################################################################
#     MLwiN MCMC Manual
#
# 3   Variance Components Models . . . . . . . . . . . . . . . . . . . . .35
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

# 3.1 A 2 level variance components model for the Tutorial dataset . . . .36

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

## The highest level comes first, then the second highest and so on
(mymodel1 <- runMLwiN(normexam ~ 1 + standlrt + (1 | school) + (1 | student), estoptions = list(EstM = 1), data = tutorial))

estimates <- mymodel1@chains
par(mfrow = c(3, 2))
plot(4501:niter(estimates), estimates[4501:niter(estimates), "deviance"], xlab = "iteration", ylab = expression(paste("Est. of deviance")), 
  type = "l")
plot(4501:niter(estimates), estimates[4501:niter(estimates), "FP_Intercept"], xlab = "iteration", ylab = expression(paste("Est. of ", 
  beta[0])), type = "l")
plot(4501:niter(estimates), estimates[4501:niter(estimates), "FP_standlrt"], xlab = "iteration", ylab = expression(paste("Est. of ", 
  beta[1])), type = "l")
plot(4501:niter(estimates), estimates[4501:niter(estimates), "RP2_var_Intercept"], xlab = "iteration", ylab = expression(paste("Est. of ", 
  sigma[u0]^2)), type = "l")
plot(4501:niter(estimates), estimates[4501:niter(estimates), "RP1_var_Intercept"], xlab = "iteration", ylab = expression(paste("Est. of ", 
  sigma[e0]^2)), type = "l")

sixway(mymodel1@chains[, "FP_standlrt", drop = FALSE], "beta_1")
sixway(mymodel1@chains[, "RP2_var_Intercept", drop = FALSE], "sigma2u0")

# 3.2 DIC and multilevel models . . . . . . . . . . . . . . . . . . . . . 41

# 3.3 Comparison between fixed and random school effects . . . . . . . . .41

(mymodel2 <- runMLwiN(normexam ~ 1 + standlrt + sex + (1 | school) + (1 | student), estoptions = list(EstM = 1), data = tutorial))

# Chapter learning outcomes . . . . . . . . . . . . . . . . . . . . . . . 43





############################################################################
