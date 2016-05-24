############################################################################
#     MLwiN User Manual
#
# 9   Logistic Models for Binary and Binomial Responses . . . . . . . . .117
#
#     Rasbash, J., Steele, F., Browne, W. J. and Goldstein, H. (2012).
#     A User's Guide to MLwiN, v2.26. Centre for Multilevel Modelling,
#     University of Bristol.
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


# 9.1 Introduction and description of the example data . . . . . . . . . 117

data(bang, package = "R2MLwiN")
summary(bang)

# 9.2 Single-level logistic regression . . . . . . . . . . . . . . . . . 119

# Link functions . . . . . . . . . . . . . . . . . . . . . . . . . . . . 119

# Interpretation of coeficients . . . . . . . . . . . . . . . . . . . . .120

# Fitting a single-level logit model in MLwiN . . . . . . . . . . . . . .121

addmargins(with(bang, table(lc, use)))

(mymodel1 <- runMLwiN(logit(use) ~ 1 + lc, D = "Binomial", data = bang))

if (!require(car)) install.packages("car")
library(car)

linearHypothesis(mymodel1, "FP_lcOne_child = FP_lcTwo_children")

# A probit model . . . . . . . . . . . . . . . . . . . . . . . . . . . . 126

(mymodel2 <- runMLwiN(probit(use) ~ 1 + lc, D = "Binomial", data = bang))

(mymodel3 <- runMLwiN(logit(use) ~ 1 + lc + age, D = "Binomial", data = bang))

# 9.3 A two-level random intercept model . . . . . . . . . . . . . . . . 128

# Model specification . . . . . . . . . . . . . . . . . . . . . . . . . .128

# Estimation procedures . . . . . . . . . . . . . . . . . . . . . . . . .128

# Fitting a two-level random intercept model in MLwiN . . . . . . . . . .129

(mymodel4 <- runMLwiN(logit(use) ~ 1 + lc + age + (1 | district), D = "Binomial", data = bang))

(mymodel5 <- runMLwiN(logit(use) ~ 1 + lc + age + (1 | district), D = "Binomial", estoptions = list(nonlinear = c(N = 1,
  M = 2), startval = list(FP.b = mymodel4@FP, FP.v = mymodel4@FP.cov, RP.b = mymodel4@RP, RP.v = mymodel4@RP.cov)),
  data = bang))

if (!require(car)) install.packages("car")
library(car)

linearHypothesis(mymodel5, "RP2_var_Intercept = 0")

# Variance partition coeficient . . . . . . . . . . . . . . . . . . . . .131

set.seed(1)

invlogit <- function(x) exp(x)/(1 + exp(x))

u <- sqrt(coef(mymodel5)["RP2_var_Intercept"]) * qnorm(runif(5000))

p1 <- invlogit(coef(mymodel5)["FP_Intercept"] + u)

p2 <- invlogit(coef(mymodel5)["FP_Intercept"] + coef(mymodel5)["FP_lc3plus"] + coef(mymodel5)["FP_age"] * -9.7 + u)

p3 <- invlogit(coef(mymodel5)["FP_Intercept"] + coef(mymodel5)["FP_age"] * 15.3 + u)

v1 <- p1 * (1 - p1)
lev2var1 <- sd(p1)^2
lev1var1 <- mean(v1)

v2 <- p2 * (1 - p2)
lev2var2 <- sd(p2)^2
lev1var2 <- mean(v2)

v3 <- p3 * (1 - p3)
lev2var3 <- sd(p3)^2
lev1var3 <- mean(v3)

cat(paste0("VPC = ", lev2var1/(lev2var1 + lev1var1)))

cat(paste0("VPC for a young women with 3+ children (low probability use) = ", lev2var2/(lev2var2 + lev1var2)))

cat(paste0("VPC for an old woman with no children (high probability use) = ", lev2var3/(lev2var3 + lev1var3)))

# Adding further explanatory variables . . . . . . . . . . . . . . . . . 134

table(bang$educ)

(mymodel6 <- runMLwiN(logit(use) ~ 1 + lc + age + urban + educ + hindu + (1 | district), D = "Binomial", estoptions = list(nonlinear = c(N = 1,
  M = 2), startval = list(FP.b = mymodel5@FP, FP.v = mymodel5@FP.cov, RP.b = mymodel5@RP, RP.v = mymodel5@RP.cov)),
  data = bang))

# 9.4 A two-level random coeficient model . . . . . . . . . . . . . . . .135

(mymodel7 <- runMLwiN(logit(use) ~ 1 + lc + age + urban + educ + hindu + (1 + urban | district), D = "Binomial",
  estoptions = list(nonlinear = c(N = 1, M = 2), startval = list(FP.b = mymodel6@FP, FP.v = mymodel6@FP.cov, RP.b = mymodel6@RP,
    RP.v = mymodel6@RP.cov)), data = bang))

if (!require(car)) install.packages("car")
library(car)

linearHypothesis(mymodel7, "RP2_cov_Intercept_urbanUrban = 0")
linearHypothesis(mymodel7, "RP2_var_urbanUrban = 0")
linearHypothesis(mymodel7, c("RP2_cov_Intercept_urbanUrban = 0", "RP2_var_urbanUrban = 0"))

(mymodel8 <- runMLwiN(logit(use) ~ 1 + lc + age + urban + educ + hindu + d_lit + d_pray + (1 + urban | district),
  D = "Binomial", estoptions = list(nonlinear = c(N = 1, M = 2), startval = list(FP.b = mymodel7@FP, FP.v = mymodel7@FP.cov,
    RP.b = mymodel7@RP, RP.v = mymodel7@RP.cov)), data = bang))

# 9.5 Modelling binomial data . . . . . . . . . . . . . . . . . . . . . .139

# Modelling district-level variation with district-level proportions . . 139

# Creating a district-level data set . . . . . . . . . . . . . . . . . . 140

if (!require(doBy)) install.packages("doBy")
library(doBy)

bangshort <- summaryBy(use + cons ~ district + d_lit + d_pray, FUN = c(mean, sum), data = bang)
bangshort$use.sum <- NULL
colnames(bangshort) <- c("district", "d_lit", "d_pray", "use", "cons", "denom")
bangshort$use <- bangshort$use - 1

# Fitting the model . . . . . . . . . . . . . . . . . . . . . . . . . . .142

(mymodel9 <- runMLwiN(logit(use, denom) ~ 1 + d_lit + d_pray + (1 | district), D = "Binomial", data = bangshort))

(mymodel10 <- runMLwiN(logit(use, denom) ~ 1 + d_lit + d_pray + (1 | district), D = "Binomial", estoptions = list(nonlinear = c(N = 1,
  M = 2), startval = list(FP.b = mymodel9@FP, FP.v = mymodel9@FP.cov, RP.b = mymodel9@RP, RP.v = mymodel9@RP.cov)),
  data = bangshort))

############################################################################
