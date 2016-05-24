############################################################################
#     MLwiN MCMC Manual
#
# 11  Poisson Response Modelling . . . . . . . . . . . . . . . . . . . . 153
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

# User's input if necessary

## Read mmmec data
data(mmmec, package = "R2MLwiN")

# 11.1 Simple Poisson regression model . . . . . . . . . . . . . . . . . 155

(mymodel1 <- runMLwiN(log(obs) ~ 1 + uvbi + offset(log(exp)), D = "Poisson", estoptions = list(EstM = 1, mcmcMeth = list(iterations = 50000)), 
  data = mmmec))
summary(mymodel1@chains[, "FP_uvbi"])
sixway(mymodel1@chains[, "FP_uvbi", drop = FALSE], "beta_1")

# 11.2 Adding in region level random effects . . . . . . . . . . . . . . 157

(mymodel2 <- runMLwiN(log(obs) ~ 1 + uvbi + offset(log(exp)) + (1 | region), D = "Poisson", estoptions = list(EstM = 1, 
  mcmcMeth = list(iterations = 50000, seed = 13)), data = mmmec))
summary(mymodel2@chains[, "FP_uvbi"])
sixway(mymodel2@chains[, "FP_uvbi", drop = FALSE], "beta_1")

# 11.3 Including nation effects in the model . . . . . . . . . . . . . . 159

(mymodel3 <- runMLwiN(log(obs) ~ 1 + uvbi + offset(log(exp)) + (1 | nation) + (1 | region), D = "Poisson", estoptions = list(EstM = 1, 
  mcmcMeth = list(iterations = 50000, seed = 13)), data = mmmec))


contrasts(mmmec$nation, 9) <- diag(9)

(mymodel4 <- runMLwiN(log(obs) ~ 0 + uvbi + nation + offset(log(exp)) + (1 | region), D = "Poisson", estoptions = list(EstM = 1, 
  mcmcMeth = list(iterations = 50000)), data = mmmec))

# 11.4 Interaction with UV exposure . . . . . . . . . . . . . . . . . . .161

(mymodel5 <- runMLwiN(log(obs) ~ 0 + nation + nation:uvbi + offset(log(exp)) + (1 | region), D = "Poisson", estoptions = list(EstM = 1, 
  mcmcMeth = list(iterations = 50000)), data = mmmec))
sixway(mymodel5@chains[, "FP_nationBelgium", drop = FALSE], acf.maxlag = 5000, "beta_1")

# 11.5 Problems with univariate updating Metropolis procedures . . . . . 163

(mymodel6 <- runMLwiN(log(obs) ~ 0 + nation + nation:uvbi + offset(log(exp)) + (1 | region), D = "Poisson", estoptions = list(EstM = 1, 
  mcmcMeth = list(iterations = 500000, thinning = 10)), data = mmmec))
sixway(mymodel6@chains[, "FP_nationBelgium", drop = FALSE], "beta_1")

# Chapter learning outcomes . . . . . . . . . . . . . . . . . . . . . . .128





############################################################################
