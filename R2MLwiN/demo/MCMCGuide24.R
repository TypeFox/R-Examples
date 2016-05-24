############################################################################
#     MLwiN MCMC Manual
#
# 24  Parameter expansion . . . . . . . . . . . . . . . . . . . . . . . .381
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

# 24.1 What is Parameter Expansion? . . . . . . . . . . . . . . . . . . .381

# 24.2 The tutorial example . . . . . . . . . . . . . . . . . . . . . . .383

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

## openbugs executable
if (!exists("openbugs")) openbugs <- "C:/Program Files (x86)/OpenBUGS/OpenBUGS323/OpenBUGS.exe"
while (!file.access(openbugs, mode = 0) == 0 || !file.access(openbugs, mode = 1) == 0 || !file.access(openbugs, mode = 4) == 
  0) {
  cat("Please specify the path for the OpenBUGS executable:\n")
  openbugs <- scan(what = character(0), sep = "\n")
  openbugs <- gsub("\\", "/", openbugs, fixed = TRUE)
}

## winbugs executable winbugs='C:/Program Files (x86)/WinBUGS14/WinBUGS14.exe'

## Read tutorial data
data(tutorial, package = "R2MLwiN")

## Define the model

(mymodel <- runMLwiN(normexam ~ 1 + standlrt + (1 | school) + (1 | student), estoptions = list(EstM = 1), data = tutorial))

summary(mymodel@chains[, "RP2_var_Intercept"])
sixway(mymodel@chains[, "RP2_var_Intercept", drop = FALSE], "sigma2u2")

## Parameter expansion at level 2

(mymodel <- runMLwiN(normexam ~ 1 + standlrt + (1 | school) + (1 | student), estoptions = list(EstM = 1, mcmcOptions = list(paex = c(2, 
  1))), data = tutorial))

sixway(mymodel@chains[, "RP2_var_Intercept", drop = FALSE], "sigma2u0")

# 24.3 Binary responses - Voting example . . . . . . . . . . . . . . . . 386


## Read bes83 data
data(bes83, package = "R2MLwiN")

## Define the model

(mymodel <- runMLwiN(logit(votecons) ~ 1 + defence + unemp + taxes + privat + (1 | area), D = "Binomial", estoptions = list(EstM = 1), 
  data = bes83))

sixway(mymodel@chains[, "RP2_var_Intercept", drop = FALSE], acf.maxlag = 500, "sigma2u0")

## Parameter expansion at level 2

(mymodel <- runMLwiN(logit(votecons) ~ 1 + defence + unemp + taxes + privat + (1 | area), D = "Binomial", estoptions = list(EstM = 1, 
  mcmcOptions = list(paex = c(2, 1))), data = bes83))

sixway(mymodel@chains[, "RP2_var_Intercept", drop = FALSE], acf.maxlag = 500, "sigma2u0")

# 24.4 The choice of prior distribution . . . . . . . . . . . . . . . . .390

## Uniform on the variance scale priors+Parameter expansion at level 2

(mymodel <- runMLwiN(logit(votecons) ~ 1 + defence + unemp + taxes + privat + (1 | area), D = "Binomial", estoptions = list(EstM = 1, 
  mcmcMeth = list(priorcode = 0), mcmcOptions = list(paex = c(2, 1))), data = bes83))

sixway(mymodel@chains[, "RP2_var_Intercept", drop = FALSE], acf.maxlag = 100, "sigma2u0")

# 24.5 Parameter expansion and WinBUGS . . . . . . . . . . . . . . . . . 391

mymodel <- runMLwiN(logit(votecons) ~ 1 + defence + unemp + taxes + privat + (1 | area), D = "Binomial", estoptions = list(EstM = 1, 
  mcmcMeth = list(priorcode = 0), mcmcOptions = list(paex = c(2, 1)), show.file = TRUE), BUGO = c(version = 4, n.chains = 1, 
  debug = FALSE, seed = 1, bugs = openbugs, OpenBugs = TRUE), data = bes83)

summary(mymodel)
effectiveSize(mymodel)
sixway(mymodel[, "sigma2.u2", drop = FALSE], acf.maxlag = 250)
sixway(mymodel[, "sigma2.v2", drop = FALSE], acf.maxlag = 100)

# 24.6 Parameter expansion and random slopes . . . . . . . . . . . . . . 396

## Read tutorial data
data(tutorial, package = "R2MLwiN")

## Define the model

(mymodel <- runMLwiN(normexam ~ 1 + standlrt + (1 + standlrt | school) + (1 | student), estoptions = list(EstM = 1), 
  data = tutorial))

## Parameter expansion at level 2
(mymodel <- runMLwiN(normexam ~ 1 + standlrt + (1 + standlrt | school) + (1 | student), estoptions = list(EstM = 1, 
  mcmcOptions = list(paex = c(2, 1))), data = tutorial))

# Chapter learning outcomes . . . . . . . . . . . . . . . . . . . . . . .399





############################################################################
