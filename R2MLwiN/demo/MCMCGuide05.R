############################################################################
#     MLwiN MCMC Manual
#
# 5   Prior Distributions, Starting Values and Random Number Seeds . . . .61
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

# 5.1 Prior distributions . . . . . . . . . . . . . . . . . . . . . . . . 61

# 5.2 Uniform on variance scale priors . . . . . . . . . . . . . . . . . .61

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

## IGLS
(mymodel1 <- runMLwiN(normexam ~ 1 + standlrt + (1 | school) + (1 | student), data = tutorial))

## Diffuse priors (Gamma priors)
(mymodel2 <- runMLwiN(normexam ~ 1 + standlrt + (1 | school) + (1 | student), estoptions = list(EstM = 1), data = tutorial))

## Diffuse priors (Uniform priors)
(mymodel3 <- runMLwiN(normexam ~ 1 + standlrt + (1 | school) + (1 | student), estoptions = list(EstM = 1, mcmcMeth = list(priorcode = 0)), 
  data = tutorial))

aa <- cbind(mymodel1@FP, mymodel2@FP, mymodel3@FP)
bb <- cbind(mymodel1@RP, mymodel2@RP, mymodel3@RP)
ctable <- round(rbind(aa, bb), 3)
colnames(ctable) <- c("IGLS", "Gibbs1", "Gibbs2")
print(ctable)
rm(list = c("mymodel1", "mymodel2", "mymodel3"))

# 5.3 Using informative priors . . . . . . . . . . . . . . . . . . . . . .62

## Informative normal prior for beta_1 Fit the model
(mymodel4 <- runMLwiN(normexam ~ 1 + standlrt + (1 | school) + (1 | student), estoptions = list(EstM = 1, mcmcMeth = list(priorParam = list(fixe = list(standlrt = c(1, 
  0.01))))), data = tutorial))

sixway(mymodel4@chains[, "FP_standlrt", drop = FALSE], "beta_1")

## Informative normal prior for beta_1
(mymodel5 <- runMLwiN(normexam ~ 1 + standlrt + (1 | school) + (1 | student), estoptions = list(EstM = 1, mcmcMeth = list(priorParam = list(fixe = list(standlrt = c(1, 
  0.1))))), data = tutorial))

sixway(mymodel5@chains[, "FP_standlrt", drop = FALSE], "beta_1")

# 5.4 Specifying an informative prior for a random parameter . . . . . . .65

## Specifies an informative prior for sigma2u
(mymodel6 <- runMLwiN(normexam ~ 1 + standlrt + (1 | school) + (1 | student), estoptions = list(EstM = 1, mcmcMeth = list(priorParam = list(rp2 = list(estimate = 0.2, 
  size = 100)))), data = tutorial))

sixway(mymodel6@chains[, "RP2_var_Intercept", drop = FALSE], "sigma^2_u0")

# 5.5 Changing the random number seed and the parameter starting values  .66

## Set starting values for random and fixed parameter estimates
FP.b <- c(-2, 5)
names(FP.b) <- c("FP_Intercept", "FP_standlrt")
RP.b <- c(2, 4)
names(RP.b) <- c("RP2_var_Intercept", "RP1_var_Intercept")
startval <- list(FP.b = FP.b, RP.b = RP.b)

(mymodel7 <- runMLwiN(normexam ~ 1 + standlrt + (1 | school) + (1 | student), estoptions = list(EstM = 1, mcmcMeth = list(burnin = 0, 
  iterations = 500), startval = startval), data = tutorial))

rm(list = c("mymodel4", "mymodel5", "mymodel6", "mymodel7"))

## Use different seeds
(mymodel8 <- runMLwiN(normexam ~ 1 + standlrt + (1 | school) + (1 | student), estoptions = list(EstM = 1, mcmcMeth = list(nchains=4, seed = 1:4)),
  data = tutorial))

ctable <- round(as.data.frame(lapply(mymodel8@chains[, c(names(mymodel8@FP), names(mymodel8@RP))], colMeans)), 3)
colnames(ctable) <- c("Seed1", "Seed2", "Seed3", "Seed4")
print(ctable)
rm(list = "mymodel8")

# 5.6 Improving the speed of MCMC Estimation . . . . . . . . . . . . . . .69

# Chapter learning outcomes . . . . . . . . . . . . . . . . . . . . . . . 70





############################################################################
