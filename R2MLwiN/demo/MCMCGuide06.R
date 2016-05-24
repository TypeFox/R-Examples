############################################################################
#     MLwiN MCMC Manual
#
# 6   Random Slopes Regression Models . . . . . . . . . . . . . . . . . . 71
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

## Choose MCMC algoritm for estimation (IGLS will be used to obtain starting values for MCMC)
(mymodel <- runMLwiN(normexam ~ 1 + standlrt + school + school:standlrt + (1 | student), estoptions = list(EstM = 1), 
  data = tutorial))

## Define the model Choose IGLS algoritm for estimation Fit the model
(mymodel0a <- runMLwiN(normexam ~ 1 + standlrt + (1 + standlrt | school) + (1 | student), data = tutorial))

## Choose MCMC algoritm for estimation (IGLS will be used to obtain starting values for MCMC)
(mymodel0 <- runMLwiN(normexam ~ 1 + standlrt + (1 + standlrt | school) + (1 | student), estoptions = list(EstM = 1), 
  data = tutorial))

# 6.1 Prediction intervals for a random slopes regression model . . . . . 75

## Save level 2 residual chains
(mymodel <- runMLwiN(normexam ~ 1 + standlrt + (1 + standlrt | school) + (1 | student), estoptions = list(EstM = 1, 
  mcmcMeth = list(iterations = 5001), resi.store.levs = 2), data = tutorial))

predLines(mymodel, xname = "standlrt", lev = 2, selected = NULL, probs = c(0.025, 0.975), legend.space = "right", 
  legend.ncol = 2)
dev.new()
predLines(mymodel, xname = "standlrt", lev = 2, selected = c(30, 44, 53, 59), probs = c(0.025, 0.975))

# 6.2 Alternative priors for variance matrices . . . . . . . . . . . . . .78

# 6.3 WinBUGS priors (Prior 2) . . . . . . . . . . . . . . . . . . . . . .78

## Change the starting values for Level 2 variance matrix to .1 on diagonal 0 otherwise.
RP.b <- c(0.1, 0, 0.1, 0.554)
names(RP.b) <- c("RP2_var_Intercept", "RP2_cov_Intercept_standlrt", "RP2_var_standlrt", "RP1_var_Intercept")

## Fit the model
(mymodel1 <- runMLwiN(normexam ~ 1 + standlrt + (1 + standlrt | school) + (1 | student), estoptions = list(EstM = 1, 
  startval = list(RP.b = RP.b)), data = tutorial))

# 6.4 Uniform prior . . . . . . . . . . . . . . . . . . . . . . . . . . . 79

## Diffuse priors (Uniform priors)
(mymodel2 <- runMLwiN(normexam ~ 1 + standlrt + (1 + standlrt | school) + (1 | student), estoptions = list(EstM = 1, 
  mcmcMeth = list(priorcode = 0)), data = tutorial))

# 6.5 Informative prior . . . . . . . . . . . . . . . . . . . . . . . . . 80

## Informative normal prior for Sigma_u
(mymodel3 <- runMLwiN(normexam ~ 1 + standlrt + (1 + standlrt | school) + (1 | student), estoptions = list(EstM = 1, 
  mcmcMeth = list(priorParam = list(rp2 = list(estimate = matrix(c(0.09, 0.018, 0.018, 0.015), 2, 2), size = 65)))), 
  data = tutorial))

# 6.6 Results . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 81

cat("The mean parameter estimates\n")
aa <- cbind(mymodel0a@FP, mymodel0@FP, mymodel1@FP, mymodel2@FP, mymodel3@FP)
bb <- cbind(mymodel0a@RP, mymodel0@RP, mymodel1@RP, mymodel2@RP, mymodel3@RP)
ctable <- round(rbind(aa, bb), 3)
colnames(ctable) <- c("IGLS", "default", "prior 2", "uniform", "prior 4")
print(ctable)

cat("The standard errors of parameter estimates\n")
cc <- cbind(sqrt(diag(mymodel0a@FP.cov)), sqrt(diag(mymodel0@FP.cov)), sqrt(diag(mymodel1@FP.cov)), sqrt(diag(mymodel2@FP.cov)), 
  sqrt(diag(mymodel3@FP.cov)))
dd <- cbind(sqrt(diag(mymodel0a@RP.cov)), sqrt(diag(mymodel0@RP.cov)), sqrt(diag(mymodel1@RP.cov)), sqrt(diag(mymodel2@RP.cov)), 
  sqrt(diag(mymodel3@RP.cov)))
sdtable <- round(rbind(cc, dd), 3)
colnames(sdtable) <- c("IGLS", "default", "prior 2", "uniform", "prior 4")
print(sdtable)

# Chapter learning outcomes . . . . . . . . . . . . . . . . . . . . . . . 81





############################################################################
