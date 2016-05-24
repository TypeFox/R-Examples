############################################################################
#     MLwiN MCMC Manual
#
# 23  Using Orthogonal fixed effect vectors . . . . . . . . . . . . . . .357
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

# 23.1 A simple example . . . . . . . . . . . . . . . . . . . . . . . . .358

# 23.2 Constructing orthogonal vectors . . . . . . . . . . . . . . . . . 359

# 23.3 A Binomial response example . . . . . . . . . . . . . . . . . . . 360

library(R2MLwiN)
# MLwiN folder
mlwin <- getOption("MLwiN_path")
while (!file.access(mlwin, mode = 1) == 0) {
  cat("Please specify the root MLwiN folder or the full path to the MLwiN executable:\n")
  mlwin <- scan(what = character(0), sep = "\n")
  mlwin <- gsub("\\", "/", mlwin, fixed = TRUE)
}
options(MLwiN_path = mlwin)

## openbugs executable
if (!exists("openbugs")) openbugs <- "C:/Program Files (x86)/OpenBUGS/OpenBUGS323/OpenBUGS.exe"
while (!file.access(openbugs, mode = 0) == 0 || !file.access(openbugs, mode = 1) == 0 || !file.access(openbugs, mode = 4) ==
  0) {
  cat("Please specify the path for the OpenBUGS executable:\n")
  openbugs <- scan(what = character(0), sep = "\n")
  openbugs <- gsub("\\", "/", openbugs, fixed = TRUE)
}

## winbugs executable
#winbugs="C:/Program Files (x86)/WinBUGS14/WinBUGS14.exe"

data(bang1, package="R2MLwiN")

## Define the model

(mymodel <- runMLwiN(logit(use) ~ 1 + age + lc + urban + (1 + urban | district), D = "Binomial", estoptions = list(EstM = 1),
  data = bang1))

trajectories(mymodel)

## Orthogonal update

(mymodel <- runMLwiN(logit(use) ~ 1 + age + lc + urban + (1 + urban | district), D = "Binomial", estoptions = list(EstM = 1,
  mcmcOptions = list(orth = 1)), data = bang1))

trajectories(mymodel)

# 23.4 A Poisson example . . . . . . . . . . . . . . . . . . . . . . . . 364

## Read mmmec data
data(mmmec, package = "R2MLwiN")

contrasts(mmmec$nation, 9) <- diag(9)

(mymodel <- runMLwiN(log(obs) ~ 0 + nation + nation:uvbi + offset(log(exp)) + (1 | region), D = "Poisson", estoptions = list(EstM = 1,
  mcmcMeth = list(iterations = 50000)), data = mmmec))

sixway(mymodel@chains[, "FP_nationBelgium", drop = FALSE], acf.maxlag = 5000, "beta_1")

## Orthogonal update

(mymodel <- runMLwiN(log(obs) ~ 0 + nation + nation:uvbi + offset(log(exp)) + (1 | region), D = "Poisson", estoptions = list(EstM = 1,
  mcmcMeth = list(iterations = 50000), mcmcOptions = list(orth = 1)), data = mmmec))

sixway(mymodel@chains[, "FP_nationBelgium", drop = FALSE], acf.maxlag = 100, "beta_1")

# 23.5 An Ordered multinomial example . . . . . . . . . . . . . . . . . .368

## Read alevchem data
data(alevchem, package = "R2MLwiN")

# Note: Establishment codes on their own do not uniquely identify schools.
# Schools are instead uniquely identified by LEA code, establishment ID
# combination. Thus, here we generated a unique school ID.
alevchem$school <- as.numeric(factor(paste0(alevchem$lea, alevchem$estab)))

alevchem$gcseav <- double2singlePrecision(alevchem$gcse_tot/alevchem$gcse_no - 6)

## MCMC
(mymodel <- runMLwiN(logit(a_point, cons, 6) ~ 1 + gcseav[1:5] + I(gcseav^2)[1:5] + gender[1:5] + (1[1:5] | school),
  D = "Ordered Multinomial", estoptions = list(EstM = 1), data = alevchem))

trajectories(mymodel)

## Orthogonal update
(mymodel <- runMLwiN(logit(a_point, cons, 6) ~ 1 + gcseav[1:5] + I(gcseav^2)[1:5] + gender[1:5] + (1[1:5] | school),
  D = "Ordered Multinomial", estoptions = list(EstM = 1, mcmcOptions = list(orth = 1)), data = alevchem))

trajectories(mymodel)

# 23.6 The WinBUGS interface . . . . . . . . . . . . . . . . . . . . . . 372

## Read bang1 data
data(bang1, package = "R2MLwiN")

## Orthogonal update (WinBUGS)

mymodel <- runMLwiN(logit(use) ~ 1 + age + lc + urban + (1 + urban | district), D = "Binomial", estoptions = list(EstM = 1,
  mcmcOptions = list(orth = 1), show.file = TRUE), BUGO = c(version = 4, n.chains = 1, debug = FALSE, seed = 1,
  bugs = openbugs, OpenBugs = TRUE), data = bang1)

summary(mymodel)
effectiveSize(mymodel)
sixway(mymodel[, "beta[1]", drop = FALSE])

# Chapter learning outcomes . . . . . . . . . . . . . . . . . . . . . . .379





############################################################################
