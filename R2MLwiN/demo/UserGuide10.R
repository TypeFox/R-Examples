############################################################################
#     MLwiN User Manual
#
# 10  Multinomial Logistic Models for Unordered Categorical Responses . .145
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


# 10.1 Introduction . . . . . . . . . . . . . . . . . . . . . . . . . . .145

data(bang, package = "R2MLwiN")

addmargins(with(bang, table(use4)))

# 10.2 Single-level multinomial logistic regression . . . . . . . . . . .146


# 10.3 Fitting a single-level multinomial logistic model in MLwiN . . . .147

addmargins(with(bang, table(lc, use4)))

bang$use4 <- relevel(bang$use4, 4)

(mymodel1 <- (runMLwiN(logit(use4) ~ 1 + lc, D = "Unordered Multinomial", data = bang)))

cat(paste("Pr(y = 1) =", round(exp(mymodel1@FP["FP_Intercept_Sterilization"])/(1 + exp(mymodel1@FP["FP_Intercept_Sterilization"]) + 
  exp(mymodel1@FP["FP_Intercept_Modern_reversible_method"]) + exp(mymodel1@FP["FP_Intercept_Traditional_method"])), 4), "\n"))
cat(paste("Pr(y = 2) =", round(exp(mymodel1@FP["FP_Intercept_Modern_reversible_method"])/(1 + exp(mymodel1@FP["FP_Intercept_Sterilization"]) + exp(mymodel1@FP["FP_Intercept_Modern_reversible_method"]) + 
  exp(mymodel1@FP["FP_Intercept_Traditional_method"])), 4), "\n"))
cat(paste("Pr(y = 3) =", round(exp(mymodel1@FP["FP_Intercept_Traditional_method"])/(1 + exp(mymodel1@FP["FP_Intercept_Sterilization"]) + 
  exp(mymodel1@FP["FP_Intercept_Modern_reversible_method"]) + exp(mymodel1@FP["FP_Intercept_Traditional_method"])), 4), "\n"))
cat(paste("Pr(y = 4) =", round(1/(1 + exp(mymodel1@FP["FP_Intercept_Sterilization"]) + exp(mymodel1@FP["FP_Intercept_Modern_reversible_method"]) + 
  exp(mymodel1@FP["FP_Intercept_Traditional_method"])), 4), "\n"))

# 10.4 A two-level random intercept multinomial logistic regression model 154

# 10.5 Fitting a two-level random intercept model . . . . . . . . . . . .155

(mymodel2 <- (runMLwiN(logit(use4) ~ 1 + lc + (1 | district), D = "Unordered Multinomial", data = bang)))

(mymodel3 <- (runMLwiN(logit(use4) ~ 1 + lc + (1 | district), D = "Unordered Multinomial", estoptions = list(nonlinear = c(1, 
  2), startval = list(FP.b = mymodel2@FP, FP.v = mymodel2@FP.cov, RP.b = mymodel2@RP, RP.v = mymodel2@RP.cov), resi.store = TRUE), 
  data = bang)))


mymodel3@RP["RP2_cov_Intercept_Sterilization_Intercept_Modern_reversible_method"]/sqrt(mymodel3@RP["RP2_var_Intercept_Sterilization"] * mymodel3@RP["RP2_var_Intercept_Modern_reversible_method"])
mymodel3@RP["RP2_cov_Intercept_Sterilization_Intercept_Traditional_method"]/sqrt(mymodel3@RP["RP2_var_Intercept_Sterilization"] * mymodel3@RP["RP2_var_Intercept_Traditional_method"])
mymodel3@RP["RP2_cov_Intercept_Modern_reversible_method_Intercept_Traditional_method"]/sqrt(mymodel3@RP["RP2_var_Intercept_Modern_reversible_method"] * mymodel3@RP["RP2_var_Intercept_Traditional_method"])

hipos <- rep(0, 2)
hipos[1] <- which(levels(as.factor(bang$district)) == 56)
hipos[2] <- which(levels(as.factor(bang$district)) == 11)

u0 <- mymodel3@residual$lev_2_resi_est_Intercept.Sterilization
u0se <- sqrt(mymodel3@residual$lev_2_resi_var_Intercept.Sterilization)
u0rank <- rank(u0)
u0rankhi <- u0 + u0se
u0ranklo <- u0 - u0se
u0rankno <- order(u0rank)
plot(1:60, u0[u0rankno], ylim = c(-2, 2), pch = 15, xlab = "Rank", ylab = "u0 residual estimate")
points(1:60, u0rankhi[u0rankno], pch = 24, bg = "grey")
points(1:60, u0ranklo[u0rankno], pch = 25, bg = "grey")
for (i in 1:60) lines(rep(i, 2), c(u0ranklo[u0rankno[i]], u0rankhi[u0rankno[i]]))
for (i in 1:2) points(x = which(u0rankno == hipos[i]), y = u0[u0rankno[which(u0rankno == hipos[i])]], pch = 22, bg = i + 
  1)

u1 <- mymodel3@residual$lev_2_resi_est_Intercept.Modern_reversible_method
u1se <- sqrt(mymodel3@residual$lev_2_resi_var_Intercept.Modern_reversible_method)
u1rank <- rank(u1)
u1rankhi <- u1 + u1se
u1ranklo <- u1 - u1se
u1rankno <- order(u1rank)
plot(1:60, u1[u1rankno], ylim = c(-2, 2), pch = 15, xlab = "Rank", ylab = "u1 residual estimate")
points(1:60, u1rankhi[u1rankno], pch = 24, bg = "grey")
points(1:60, u1ranklo[u1rankno], pch = 25, bg = "grey")
for (i in 1:60) lines(rep(i, 2), c(u1ranklo[u1rankno[i]], u1rankhi[u1rankno[i]]))
for (i in 1:2) points(x = which(u1rankno == hipos[i]), y = u1[u1rankno[which(u1rankno == hipos[i])]], pch = 22, bg = i + 
  1)

u2 <- mymodel3@residual$lev_2_resi_est_Intercept.Traditional_method
u2se <- sqrt(mymodel3@residual$lev_2_resi_var_Intercept.Traditional_method)
u2rank <- rank(u2)
u2rankhi <- u2 + u2se
u2ranklo <- u2 - u2se
u2rankno <- order(u2rank)
plot(1:60, u2[u2rankno], ylim = c(-2, 2), pch = 15, xlab = "Rank", ylab = "u2 residual estimate")
points(1:60, u2rankhi[u2rankno], pch = 24, bg = "grey")
points(1:60, u2ranklo[u2rankno], pch = 25, bg = "grey")
for (i in 1:60) lines(rep(i, 2), c(u2ranklo[u2rankno[i]], u2rankhi[u2rankno[i]]))
for (i in 1:2) points(x = which(u2rankno == hipos[i]), y = u2[u2rankno[which(u2rankno == hipos[i])]], pch = 22, bg = i + 
  1)

# Chapter learning outcomes . . . . . . . . . . . . . . . . . . . . . . .159

############################################################################
