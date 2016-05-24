############################################################################
#     MLwiN MCMC Manual
#
# 19  Mixed Response Models and Correlated Residuals . . . . . . . . . . 287
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

# 19.1 Mixed response models . . . . . . . . . . . . . . . . . . . . . . 287

# 19.2 The JSP mixed response example . . . . . . . . . . . . . . . . . .289

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

## Read jspmix1 data
data(jspmix1, package = "R2MLwiN")

tab1 <- matrix(, 3, 3)
colnames(tab1) <- c("0", "1", "TOTALS")
rownames(tab1) <- c("N", "MEANS", "SDs")
tab1[1, 1:2] <- colSums(table(jspmix1$english, jspmix1$behaviour))
tab1[1, 3] <- sum(tab1[1, 1:2])
tab1[2, 1] <- mean(jspmix1$english[jspmix1$behaviour == 0])
tab1[2, 2] <- mean(jspmix1$english[jspmix1$behaviour == 1])
tab1[2, 3] <- mean(jspmix1$english)
tab1[3, 1] <- sd(jspmix1$english[jspmix1$behaviour == 0])
tab1[3, 2] <- sd(jspmix1$english[jspmix1$behaviour == 1])
tab1[3, 3] <- sd(jspmix1$english)
formatC(tab1)

vars <- cbind(as.numeric(jspmix1$sex) - 1, jspmix1$fluent, jspmix1$ravens, jspmix1$english, as.numeric(jspmix1$behaviour) - 
  1)
colnames(vars) <- c("sex", "fluent", "ravens", "english", "behaviour")
round(cor(vars), 4)

# 19.3 Setting up a single level mixed response model . . . . . . . . . .291

(mymodel <- runMLwiN(c(english, probit(behaviour)) ~ 1 + sex + ravens + fluent[1] + (1[1] | id), D = c("Mixed", 
  "Normal", "Binomial"), estoptions = list(EstM = 1, mcmcMeth = list(fixM = 1, residM = 1, Lev1VarM = 1), sort.ignore = TRUE), 
  data = jspmix1))

# 19.4 Multilevel mixed response model . . . . . . . . . . . . . . . . . 294

(mymodel <- runMLwiN(c(english, probit(behaviour)) ~ 1 + sex + ravens + fluent[1] + (1 | school) + (1[1] | id), 
  D = c("Mixed", "Normal", "Binomial"), estoptions = list(EstM = 1, mcmcMeth = list(fixM = 1, residM = 1, Lev1VarM = 1)), 
  data = jspmix1))

# 19.5 Rats dataset . . . . . . . . . . . . . . . . . . . . . . . . . . .295

## Read rats data
data(rats, package = "R2MLwiN")

(mymodel <- runMLwiN(c(y8, y15, y22, y29, y36) ~ 1 + (1 | rat), D = "Multivariate Normal", estoptions = list(EstM = 1), 
  data = rats))

sixway(mymodel@chains[, "RP1_var_Intercept_y8", drop = FALSE], "sigma2u0")

covM1 <- matrix(, 5, 5)
colnames(covM1) <- rownames(covM1) <- c("cons.y8", "cons.y15", "cons.y22", "cons.y29", "cons.y36")
covM1[upper.tri(covM1, diag = TRUE)] <- mymodel@RP
# covM1[lower.tri(covM1)] <- t(covM1)[lower.tri(covM1)]
round(t(covM1), 3)
round(cov2cor(t(covM1)), 3)

# 19.6 Fitting an autoregressive structure to the variance matrix . . . .298

(mymodel <- runMLwiN(c(y8, y15, y22, y29, y36) ~ 1 + (1 | rat), D = "Multivariate Normal", estoptions = list(EstM = 1, 
  mcmcMeth = list(iterations = 50000), mcmcOptions = list(mcco = 4)), data = rats))

covM2 <- matrix(, 5, 5)
colnames(covM2) <- rownames(covM2) <- c("cons.y8", "cons.y15", "cons.y22", "cons.y29", "cons.y36")
covM2[upper.tri(covM2, diag = TRUE)] <- mymodel@RP
# covM2[lower.tri(covM2)] <- t(covM2)[lower.tri(covM2)]
round(t(covM2), 3)
round(cov2cor(t(covM2)), 3)

# Chapter learning outcomes . . . . . . . . . . . . . . . . . . . . . . .128





############################################################################
