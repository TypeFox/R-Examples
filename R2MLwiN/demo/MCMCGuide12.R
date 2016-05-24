############################################################################
#     MLwiN MCMC Manual
#
# 12  Unordered Categorical Responses . . . . . . . . . . . . . . . . . .167
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


## Read bang data
data(bang, package = "R2MLwiN")

bang$use4 <- relevel(bang$use4, 4)

# 12.1 Fitting a first single-level multinomial model . . . . . . . . . .169

(mymodel <- runMLwiN(log(use4) ~ 1, D = "Unordered Multinomial", estoptions = list(EstM = 1), data = bang))

cat(paste("Pr(y = 1) =", round(exp(mymodel@FP["FP_Intercept_Sterilization"])/(1 + exp(mymodel@FP["FP_Intercept_Sterilization"]) + exp(mymodel@FP["FP_Intercept_Modern_reversible_method"]) + 
  exp(mymodel@FP["FP_Intercept_Traditional_method"])), 4), "\n"))
cat(paste("Pr(y = 2) =", round(exp(mymodel@FP["FP_Intercept_Modern_reversible_method"])/(1 + exp(mymodel@FP["FP_Intercept_Sterilization"]) + exp(mymodel@FP["FP_Intercept_Modern_reversible_method"]) + 
  exp(mymodel@FP["FP_Intercept_Traditional_method"])), 4), "\n"))
cat(paste("Pr(y = 3) =", round(exp(mymodel@FP["FP_Intercept_Traditional_method"])/(1 + exp(mymodel@FP["FP_Intercept_Sterilization"]) + exp(mymodel@FP["FP_Intercept_Modern_reversible_method"]) + 
  exp(mymodel@FP["FP_Intercept_Traditional_method"])), 4), "\n"))

# 12.2 Adding predictor variables . . . . . . . . . . . . . . . . . . . .173

(mymodel <- runMLwiN(log(use4) ~ 1 + lc, D = "Unordered Multinomial", estoptions = list(EstM = 1), data = bang))

cat(paste("Pr(y = 3) =", round(exp(mymodel@FP["FP_Intercept_Traditional_method"])/(1 + exp(mymodel@FP["FP_Intercept_Sterilization"]) + exp(mymodel@FP["FP_Intercept_Modern_reversible_method"]) + 
  exp(mymodel@FP["FP_Intercept_Traditional_method"])), 4), "\n"))
cat(paste("Pr(y = 3) =", round(exp(mymodel@FP["FP_Intercept_Traditional_method"] + mymodel@FP["FP_lcTwo_children_Traditional_method"])/(1 + exp(mymodel@FP["FP_Intercept_Sterilization"] + 
  mymodel@FP["FP_lcTwo_children_Sterilization"]) + exp(mymodel@FP["FP_Intercept_Modern_reversible_method"] + mymodel@FP["FP_lcTwo_children_Modern_reversible_method"]) + exp(mymodel@FP["FP_Intercept_Traditional_method"] + 
  mymodel@FP["FP_lcTwo_children_Traditional_method"])), 4), "\n"))

# 12.3 Interval estimates for conditional probabilities . . . . . . . . .175

chains <- mymodel@chains
pred1 <- exp(chains[, "FP_Intercept_Traditional_method"])/(1 + exp(chains[, "FP_Intercept_Sterilization"]) + exp(chains[, "FP_Intercept_Modern_reversible_method"]) + 
  exp(chains[, "FP_Intercept_Traditional_method"]))
summary(pred1)
sixway(pred1, "prob1")

pred2 <- exp(chains[, "FP_Intercept_Traditional_method"] + chains[, "FP_lcTwo_children_Traditional_method"])/(1 + exp(chains[, "FP_Intercept_Sterilization"] + chains[, 
  "FP_lcTwo_children_Sterilization"]) + exp(chains[, "FP_Intercept_Modern_reversible_method"] + chains[, "FP_lcTwo_children_Modern_reversible_method"]) + exp(chains[, "FP_Intercept_Traditional_method"] + 
  chains[, "FP_lcTwo_children_Traditional_method"]))
summary(pred2)
sixway(pred2, "prob1")

# 12.4 Adding district level random effects . . . . . . . . . . . . . . .177

## Uses IGLS
(mymodel <- runMLwiN(log(use4) ~ 1 + lc + (1 | district), D = "Unordered Multinomial", estoptions = list(EstM = 0, 
  nonlinear = c(1, 2)), data = bang))

## Uses MCMC
(mymodel <- runMLwiN(log(use4) ~ 1 + lc + (1 | district), D = "Unordered Multinomial", estoptions = list(EstM = 1, 
  nonlinear = c(1, 2)), data = bang))
sixway(mymodel@chains[, "RP2_var_Intercept_Sterilization", drop = FALSE], "sigma2v0")

RP3.cons <- matrix(, 3, 3)
RP3.cons[upper.tri(RP3.cons, diag = TRUE)] <- mymodel@RP[1:6]
RP3.cons[lower.tri(RP3.cons)] <- RP3.cons[upper.tri(RP3.cons)]
round(cov2cor(RP3.cons), 3)

# Chapter learning outcomes . . . . . . . . . . . . . . . . . . . . . . .128





############################################################################
