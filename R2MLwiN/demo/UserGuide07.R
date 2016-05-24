############################################################################
#     MLwiN User Manual
#
# 7   Modelling the Variance as a Function of Explanatory Variables . . . 89
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


# 7.1 A level 1 variance function for two groups . . . . . . . . . . . . .89

data(tutorial, package = "R2MLwiN")

covmatrix <- matrix(, nrow = 3, ncol = 1)
covmatrix[1, 1] <- 1
covmatrix[2, 1] <- "sexboy"
covmatrix[3, 1] <- "sexgirl"

contrasts(tutorial$sex, 2) <- diag(2)

(mymodel1 <- runMLwiN(normexam ~ 0 + sex + (0 + sex | student), estoptions = list(clre = covmatrix), data = tutorial))

contrasts(tutorial$sex, 1) <- c(0, 1)

# 7.2 Variance functions at level 2 . . . . . . . . . . . . . . . . . . . 95

(mymodel2 <- runMLwiN(normexam ~ 1 + standlrt + (1 + standlrt | school) + (1 | student), data = tutorial))

l2varfn <- mymodel2@RP["RP2_var_Intercept"] + (2 * mymodel2@RP["RP2_cov_Intercept_standlrt"] * mymodel2@data$standlrt) + 
  (mymodel2@RP["RP2_var_standlrt"] * mymodel2@data$standlrt^2)

varfndata <- as.data.frame(cbind(mymodel2@data$standlrt, l2varfn)[order(mymodel2@data$standlrt), ])
colnames(varfndata) <- c("standlrt", "l2varfn")

plot(varfndata$standlrt, varfndata$l2varfn, type = "l")

# 7.3 Further elaborating the model for the student-level variance . . . .99

(mymodel3 <- runMLwiN(normexam ~ 1 + standlrt + (1 + standlrt | school) + (1 + standlrt | student), data = tutorial))

l2varfn <- mymodel3@RP["RP2_var_Intercept"] + (2 * mymodel3@RP["RP2_cov_Intercept_standlrt"] * mymodel3@data$standlrt) + 
  (mymodel3@RP["RP2_var_standlrt"] * mymodel3@data$standlrt^2)

l1varfn <- mymodel3@RP["RP1_var_Intercept"] + (2 * mymodel3@RP["RP1_cov_Intercept_standlrt"] * mymodel3@data$standlrt) + 
  (mymodel3@RP["RP1_var_standlrt"] * mymodel3@data$standlrt^2)

varfndata <- as.data.frame(cbind(mymodel3@data$standlrt, l2varfn, l1varfn)[order(mymodel3@data$standlrt), ])
colnames(varfndata) <- c("standlrt", "l2varfn", "l1varfn")

xyplot(l2varfn + l1varfn ~ standlrt, data = varfndata, type = "l")


covmatrix <- matrix(, nrow = 3, ncol = 3)
covmatrix[1, 1] <- 1
covmatrix[2, 1] <- "standlrt"
covmatrix[3, 1] <- "standlrt"
covmatrix[1, 2] <- 1
covmatrix[2, 2] <- "sexgirl"
covmatrix[3, 2] <- "Intercept"
covmatrix[1, 3] <- 1
covmatrix[2, 3] <- "standlrt"
covmatrix[3, 3] <- "sexgirl"

(mymodel4 <- runMLwiN(normexam ~ 1 + standlrt + sex + (1 + standlrt | school) + (1 + standlrt + sex | student), estoptions = list(clre = covmatrix), 
  data = tutorial))

covmatrix <- matrix(, nrow = 3, ncol = 2)
covmatrix[1, 1] <- 1
covmatrix[2, 1] <- "standlrt"
covmatrix[3, 1] <- "standlrt"
covmatrix[1, 2] <- 1
covmatrix[2, 2] <- "sexgirl"
covmatrix[3, 2] <- "Intercept"

(mymodel5 <- runMLwiN(normexam ~ 1 + standlrt + sex + (1 + standlrt | school) + (1 + standlrt + sex | student), estoptions = list(clre = covmatrix), 
  data = tutorial))

l2varfn <- mymodel5@RP["RP2_var_Intercept"] + (2 * mymodel5@RP["RP2_cov_Intercept_standlrt"] * mymodel5@data$standlrt) + 
  (mymodel5@RP["RP2_var_standlrt"] * mymodel5@data$standlrt^2)

l1varfnboys <- mymodel5@RP["RP1_var_Intercept"] + (2 * mymodel5@RP["RP1_cov_Intercept_standlrt"] * mymodel5@data$standlrt)

l1varfngirls <- mymodel5@RP["RP1_var_Intercept"] + (2 * mymodel5@RP["RP1_cov_Intercept_standlrt"] * mymodel5@data$standlrt) + 
  (2 * mymodel5@RP["RP1_cov_standlrt_sexgirl"] * mymodel5@data$standlrt) + mymodel5@RP["RP1_var_sexgirl"]

varfndata <- as.data.frame(cbind(mymodel5@data$standlrt, l2varfn, l1varfnboys, l1varfngirls)[order(mymodel5@data$standlrt), 
  ])
colnames(varfndata) <- c("standlrt", "l2varfn", "l1varfnboys", "l1varfngirls")

xyplot(l2varfn + l1varfnboys + l1varfngirls ~ standlrt, data = varfndata, type = "l")

#     Chapter learning outcomes . . . . . . . . . . . . . . . . . . . . .106

############################################################################
