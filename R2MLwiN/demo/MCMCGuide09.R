############################################################################
#     MLwiN MCMC Manual
#
# 9   Modelling Complex Variance at Level 1 / Heteroscedasticity. . . . .111
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

## Read tutorial data
data(tutorial, package = "R2MLwiN")

boy.normexam <- tutorial$normexam[which(tutorial$sex == "boy")]
girl.normexam <- tutorial$normexam[which(tutorial$sex == "girl")]
tab1 <- cbind(c(length(boy.normexam), mean(boy.normexam), sd(boy.normexam)), c(length(girl.normexam), mean(girl.normexam), 
  sd(girl.normexam)), c(length(tutorial$normexam), mean(tutorial$normexam), sd(tutorial$normexam)))
colnames(tab1) <- c("0", "1", "TOTAL")
rownames(tab1) <- c("N", "MEANS", "SDs")
formatC(round(tab1, 6))

c5 <- tutorial$standlrt
intakecat <- rep(0, length(c5))
intakecat[which(c5 > -1)] <- 1
intakecat[which(c5 > -0.5)] <- 2
intakecat[which(c5 > -0.1)] <- 3
intakecat[which(c5 > 0.3)] <- 4
intakecat[which(c5 > 0.7)] <- 5
intakecat[which(c5 > 1.1)] <- 6
normexam <- tutorial$normexam
tab2 <- cbind(c(sum(intakecat == 0), mean(normexam[intakecat == 0]), sd(normexam[intakecat == 0])), c(sum(intakecat == 
  1), mean(normexam[intakecat == 1]), sd(normexam[intakecat == 1])), c(sum(intakecat == 2), mean(normexam[intakecat == 
  2]), sd(normexam[intakecat == 2])), c(sum(intakecat == 3), mean(normexam[intakecat == 3]), sd(normexam[intakecat == 
  3])), c(sum(intakecat == 4), mean(normexam[intakecat == 4]), sd(normexam[intakecat == 4])), c(sum(intakecat == 
  5), mean(normexam[intakecat == 5]), sd(normexam[intakecat == 5])), c(sum(intakecat == 6), mean(normexam[intakecat == 
  6]), sd(normexam[intakecat == 6])), c(length(intakecat), mean(normexam), sd(normexam)))
colnames(tab2) <- c("0", "1", "2", "3", "4", "5", "6", "TOTAL")
rownames(tab2) <- c("N", "MEANS", "SDs")
formatC(round(tab2, 6))

# 9.1 MCMC algorithm for a 1 level Normal model with complex variation . 113

# 9.2 Setting up the model in MLwiN . . . . . . . . . . . . . . . . . . .115

## The highest level comes first, then the second highest and so on
## Fit the model
(mymodel1 <- runMLwiN(normexam ~ 1 + standlrt + (1 + standlrt | student), estoptions = list(EstM = 1), data = tutorial))
trajectories(mymodel1, Range = c(4501, 5000))

l1varfn <- mymodel1@RP["RP1_var_Intercept"] + 2 * mymodel1@RP["RP1_cov_Intercept_standlrt"] * tutorial$standlrt + 
  mymodel1@RP["RP1_var_standlrt"] * tutorial$standlrt^2
plot(sort(tutorial$standlrt), l1varfn[order(tutorial$standlrt)], xlab = "standlrt", ylab = "l1varfn", type = "l")
abline(v = 0, lty = "dotted")

# 9.3 Complex variance functions in multilevel models . . . . . . . . . .119

## The highest level comes first, then the second highest and so on Choose option(s) for inference Fit the model
(mymodel2 <- runMLwiN(normexam ~ 1 + standlrt + (1 + standlrt | school) + (1 | student), estoptions = list(EstM = 1), 
  data = tutorial))

l2varfn <- mymodel2@RP["RP2_var_Intercept"] + 2 * mymodel2@RP["RP2_cov_Intercept_standlrt"] * tutorial$standlrt + 
  mymodel2@RP["RP2_var_standlrt"] * tutorial$standlrt^2
l1varfn <- mymodel2@RP["RP1_var_Intercept"]
plot(sort(tutorial$standlrt), l2varfn[order(tutorial$standlrt)], xlab = "standlrt", ylab = "varfns", ylim = c(0, 0.6), 
  type = "l")
abline(h = l1varfn)
abline(v = 0, lty = "dotted")

## Fit the model
(mymodel3 <- runMLwiN(normexam ~ 1 + standlrt + (1 + standlrt | school) + (1 + standlrt | student), estoptions = list(EstM = 1), 
  data = tutorial))

## Remove term standlrt/standlrt from the level 1 covariance matrix
clre <- matrix(, nrow = 3, ncol = 1)
clre[1, 1] <- 1
clre[2, 1] <- "standlrt"
clre[3, 1] <- "standlrt"

## Fit the model
(mymodel4 <- runMLwiN(normexam ~ 1 + standlrt + (1 + standlrt | school) + (1 + standlrt | student), estoptions = list(EstM = 1, 
  clre = clre), data = tutorial))

# 9.4 Relationship with gender . . . . . . . . . . . . . . . . . . . . . 123

tutorial$girl <- as.integer(tutorial$sex) - 1

## Remove term standlrt/standlrt and girl/girl from the level 1 covariance matrix
clre <- matrix(, nrow = 3, ncol = 2)
clre[1, 1] <- 1
clre[2, 1] <- "standlrt"
clre[3, 1] <- "standlrt"
clre[1, 2] <- 1
clre[2, 2] <- "girl"
clre[3, 2] <- "girl"

(mymodel5 <- runMLwiN(normexam ~ 1 + standlrt + girl + (1 + standlrt | school) + (1 + standlrt + girl | student), 
  estoptions = list(EstM = 1, clre = clre), data = tutorial))

l2varfn <- mymodel5@RP["RP2_var_Intercept"] + 2 * mymodel5@RP["RP2_cov_Intercept_standlrt"] * tutorial$standlrt + 
  mymodel5@RP["RP2_var_standlrt"] * tutorial$standlrt^2
l1varfnboys <- mymodel5@RP["RP1_var_Intercept"] + 2 * mymodel5@RP["RP1_cov_Intercept_standlrt"] * tutorial$standlrt
l1varfngirls <- mymodel5@RP["RP1_var_Intercept"] + 2 * mymodel5@RP["RP1_cov_Intercept_standlrt"] * tutorial$standlrt + 
  2 * mymodel5@RP["RP1_cov_Intercept_girl"] + 2 * mymodel5@RP["RP1_cov_standlrt_girl"] * tutorial$standlrt
plot(sort(tutorial$standlrt), l2varfn[order(tutorial$standlrt)], xlab = "standlrt", ylab = "varfns", ylim = c(0, 0.8), 
  type = "l")
lines(sort(tutorial$standlrt), l1varfnboys[order(tutorial$standlrt)])
lines(sort(tutorial$standlrt), l1varfngirls[order(tutorial$standlrt)])
abline(v = 0, lty = "dotted")

# 9.5 Alternative log precision formulation . . . . . . . . . . . . . . .126

(mymodel6 <- runMLwiN(normexam ~ 1 + standlrt + girl + (1 + standlrt | school) + (1 + standlrt + girl | student), 
  estoptions = list(EstM = 1, clre = clre, mcmcMeth = list(lclo = 1)), data = tutorial))

l2varfn <- mymodel6@RP["RP2_var_Intercept"] + 2 * mymodel6@RP["RP2_cov_Intercept_standlrt"] * tutorial$standlrt + 
  mymodel6@RP["RP2_var_standlrt"] * tutorial$standlrt^2
l1varfnboys <- 1/exp(mymodel6@RP["RP1_var_Intercept"] + 2 * mymodel6@RP["RP1_cov_Intercept_standlrt"] * tutorial$standlrt)
l1varfngirls <- 1/exp(mymodel6@RP["RP1_var_Intercept"] + 2 * mymodel6@RP["RP1_cov_Intercept_standlrt"] * tutorial$standlrt + 
  2 * mymodel6@RP["RP1_cov_Intercept_girl"] + 2 * mymodel6@RP["RP1_cov_standlrt_girl"] * tutorial$standlrt)
plot(sort(tutorial$standlrt), l2varfn[order(tutorial$standlrt)], xlab = "standlrt", ylab = "varfns", ylim = c(0, 0.8), 
  type = "l")
lines(sort(tutorial$standlrt), l1varfnboys[order(tutorial$standlrt)])
lines(sort(tutorial$standlrt), l1varfngirls[order(tutorial$standlrt)])
abline(v = 0, lty = "dotted")

# Chapter learning outcomes . . . . . . . . . . . . . . . . . . . . . . .128





############################################################################
