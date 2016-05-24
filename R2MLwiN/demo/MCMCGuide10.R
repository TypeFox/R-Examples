############################################################################
#     MLwiN MCMC Manual
#
# 10  Modelling Binary Responses . . . . . . . . . . . . . . . . . . . . 129
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

## Read bang1 data
data(bang1, package = "R2MLwiN")

## openbugs executable
if (!exists("openbugs")) openbugs <- "C:/Program Files (x86)/OpenBUGS/OpenBUGS323/OpenBUGS.exe"
while (!file.access(openbugs, mode = 0) == 0 || !file.access(openbugs, mode = 1) == 0 || !file.access(openbugs, mode = 4) == 
  0) {
  cat("Please specify the path for the OpenBUGS executable:\n")
  openbugs <- scan(what = character(0), sep = "\n")
  openbugs <- gsub("\\", "/", openbugs, fixed = TRUE)
}

# User's input if necessary

## winbugs executable
#winbugs="C:/Program Files (x86)/WinBUGS14/WinBUGS14.exe"

# 10.1 Simple logistic regression model . . . . . . . . . . . . . . . . .130

(mymodel1 <- runMLwiN(logit(use) ~ 1 + age, D = "Binomial", estoptions = list(EstM = 1), data = bang1))
summary(mymodel1@chains[, "FP_age"])
sixway(mymodel1@chains[, "FP_age", drop = FALSE], "beta_1")

## 15,000 iterations
(mymodel2 <- runMLwiN(logit(use) ~ 1 + age, D = "Binomial", estoptions = list(EstM = 1, mcmcMeth = list(iterations = 15000)), 
  data = bang1))
sixway(mymodel1@chains[, "FP_age", drop = FALSE], "beta_1")

## Change to 5000 iterations by default
(mymodel3 <- runMLwiN(logit(use) ~ 1 + age + lc, D = "Binomial", estoptions = list(EstM = 1), data = bang1))

# 10.2 Random effects logistic regression model . . . . . . . . . . . . .136

(mymodel4 <- runMLwiN(logit(use) ~ 1 + age + lc + (1 | district), D = "Binomial", estoptions = list(EstM = 1), 
  data = bang1))
summary(mymodel4@chains[, "RP2_var_Intercept"])
sixway(mymodel4@chains[, "RP2_var_Intercept", drop = FALSE], "sigma2u0")

# 10.3 Random coefficients for area type . . . . . . . . . . . . . . . . 139

(mymodel5 <- runMLwiN(logit(use) ~ 1 + age + lc + urban + (1 | district), D = "Binomial", estoptions = list(EstM = 1), 
  data = bang1))

(mymodel6 <- runMLwiN(logit(use) ~ 1 + age + lc + urban + (1 + urban | district), D = "Binomial", estoptions = list(EstM = 1), 
  data = bang1))

# 10.4 Probit regression . . . . . . . . . . . . . . . . . . . . . . . . 141

# 10.5 Running a probit regression in MLwiN . . . . . . . . . . . . . . .142

## Gibbs
(mymodel7 <- runMLwiN(probit(use) ~ 1 + age + lc + urban + (1 + urban | district), D = "Binomial", estoptions = list(EstM = 1, 
  mcmcMeth = list(fixM = 1, residM = 1)), data = bang1))

## Univariate MH by default
(mymodel8 <- runMLwiN(probit(use) ~ 1 + age + lc + urban + (1 + urban | district), D = "Binomial", estoptions = list(EstM = 1), 
  data = bang1))

cat("The mean parameter estimates\n")
aa <- cbind(mymodel7@FP, mymodel8@FP)
ESS.aa <- effectiveSize(mymodel7@chains[, 2:11])
bb <- cbind(mymodel7@RP, mymodel8@RP)
ESS.bb <- effectiveSize(mymodel8@chains[, 2:11])
ctable <- round(rbind(aa, bb), 3)
ctable <- cbind(ctable[, 1], round(ESS.aa), ctable[, 2], round(ESS.bb))
colnames(ctable) <- c("Gibbs", "ESS(Gibbs)", "Metropolis", "ESS(Metropolis)")
print(ctable)

cat("The standard errors of parameter estimates\n")
cc <- cbind(sqrt(diag(mymodel7@FP.cov)), sqrt(diag(mymodel8@FP.cov)))
dd <- cbind(sqrt(diag(mymodel7@RP.cov)), sqrt(diag(mymodel8@RP.cov)))
sdtable <- round(rbind(cc, dd), 3)
colnames(sdtable) <- c("Gibbs", "Metropolis")
print(sdtable)


# 10.6 Comparison with WinBUGS . . . . . . . . . . . . . . . . . . . . . 144

mymodel9 <- runMLwiN(logit(use) ~ 1 + age + (1 | district), D = "Binomial", estoptions = list(EstM = 1), BUGO = c(version = 4, 
  n.chains = 1, debug = FALSE, seed = 1, bugs = openbugs, OpenBugs = TRUE), data = bang1)

summary(mymodel9)
summary(mymodel9[, "beta[1]"])
sixway(mymodel9[, "beta[1]", drop = FALSE])

(mymodel10 <- runMLwiN(logit(use) ~ 1 + age + (1 | district), D = "Binomial", estoptions = list(EstM = 1), 
  data = bang1))

summary(mymodel10@chains[, "FP_Intercept"])
sixway(mymodel10@chains[, "FP_Intercept", drop = FALSE], "beta_0")

# Chapter learning outcomes . . . . . . . . . . . . . . . . . . . . . . .128





############################################################################  
