############################################################################
#     MLwiN User Manual
#
# 13  Fitting Models to Repeated Measures Data . . . . . . . . . . . . . 191
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


# 13.1 Introduction . . . . . . . . . . . . . . . . . . . . . . . . . . .191

# 13.2 A basic model . . . . . . . . . . . . . . . . . . . . . . . . . . 194

data(reading1, package = "R2MLwiN")
summary(reading1)

reading1[reading1 == -10] <- NA

summary(reading1)

reading <- reshape(reading1, idvar = "student", timevar = "id", varying = c("read1", "age1", "read2", "age2", "read3", 
  "age3", "read4", "age4", "read5", "age5", "read6", "age6"), sep = "", direction = "long")

reading <- reading[c("student", "id", "age", "read")]
reading <- reading[order(reading$student, reading$id), ]

colnames(reading) <- c("student", "occasion", "age", "reading")
rownames(reading) <- NULL


summary(reading)

head(reading, 5)

tab <- aggregate(reading ~ occasion, reading, function(x) c(N = length(x), mean = mean(x), sd = sd(x)))
tab <- rbind(tab, c(NA, NA))
tab$reading[7, ] <- c(length(na.omit(reading$reading)), mean(na.omit(reading$reading)), sd(na.omit(reading$reading)))
rownames(tab)[7] <- "Total"
tab

tab <- aggregate(age ~ occasion, reading, function(x) c(N = length(x), mean = mean(x), sd = sd(x)))
tab <- rbind(tab, c(NA, NA))
tab$age[7, ] <- c(length(na.omit(reading$age)), mean(na.omit(reading$age)), sd(na.omit(reading$age)))
rownames(tab)[7] <- "Total"
tab

(mymodel1 <- runMLwiN(reading ~ 1 + (1 | student) + (1 | occasion), data = reading))


# 13.3 A linear growth curve model . . . . . . . . . . . . . . . . . . . 201

(mymodel2 <- runMLwiN(reading ~ 1 + age + (1 | student) + (1 | occasion), estoptions = list(startval = list(FP.b = mymodel1@FP, 
  FP.v = mymodel1@FP.cov, RP.b = mymodel1@RP, RP.v = mymodel1@RP.cov)), data = reading))

(mymodel3 <- runMLwiN(reading ~ 1 + age + (1 + age | student) + (1 | occasion), estoptions = list(resi.store = TRUE, 
  startval = list(FP.b = mymodel2@FP, FP.v = mymodel2@FP.cov, RP.b = mymodel2@RP, RP.v = mymodel2@RP.cov)), data = reading))

u0 <- mymodel3@residual$lev_2_resi_est_Intercept
u0std <- (u0 - mean(u0))/sd(u0)

u1 <- mymodel3@residual$lev_2_resi_est_age
u1std <- (u1 - mean(u1))/sd(u1)

plot(u0std, u1std, asp = 1)

e0 <- na.omit(mymodel3@residual$lev_1_resi_est_Intercept)
e0std <- (e0 - mean(e0))/sd(e0)
e0rank <- rank(e0)
e0uniform <- (e0rank - 0.5)/length(e0rank)
e0nscore <- qnorm(e0uniform)

plot(e0std, e0nscore, asp = 1)

# 13.4 Complex level 1 variation . . . . . . . . . . . . . . . . . . . . 204

(mymodel4 <- runMLwiN(reading ~ 1 + age + (1 + age | student) + (1 + age | occasion), data = reading))

# 13.5 Repeated measures modelling of non-linear polynomial growth . . . 205

(mymodel5 <- runMLwiN(reading ~ 1 + age + I(age^2) + (1 + age + I(age^2) | student) + (1 + age | occasion), estoptions = list(resi.store = TRUE), 
  data = reading))

l2varfn <- mymodel5@RP["RP2_var_Intercept"] + 2 * mymodel5@RP["RP2_cov_Intercept_age"] * mymodel5@data$age + mymodel5@RP["RP2_var_age"] * 
  mymodel5@data$age^2 + 2 * mymodel5@RP["RP2_cov_Intercept_I(age^2)"] * mymodel5@data$age^2 + 2 * mymodel5@RP["RP2_cov_age_I(age^2)"] * 
  mymodel5@data$age * mymodel5@data[["I(age^2)"]] + mymodel5@RP["RP2_var_I(age^2)"] * mymodel5@data[["I(age^2)"]]^2

l1varfn <- mymodel5@RP["RP1_var_Intercept"] + 2 * mymodel5@RP["RP1_cov_Intercept_age"] * mymodel5@data$age + mymodel5@RP["RP1_var_age"] * 
  mymodel5@data$age^2

totvarfn <- l2varfn + l1varfn

plot(mymodel5@data$age, totvarfn)

xb <- predict(mymodel5)

u0 <- mymodel5@residual$lev_2_resi_est_Intercept
u1 <- mymodel5@residual$lev_2_resi_est_age
u2 <- mymodel5@residual[["lev_2_resi_est_I(age^2)"]]

yhat <- xb + u0[mymodel5@data$student] + u1[mymodel5@data$student] * mymodel5@data$age + u2[mymodel5@data$student] * 
  mymodel5@data[["I(age^2)"]]

plot(mymodel5@data$age, yhat, type = "n")
lines(mymodel5@data$age[mymodel5@data$student == 1], yhat[mymodel5@data$student == 1], col = 1)
lines(mymodel5@data$age[mymodel5@data$student == 2], yhat[mymodel5@data$student == 2], col = 2)
lines(mymodel5@data$age[mymodel5@data$student == 3], yhat[mymodel5@data$student == 3], col = 3)
lines(mymodel5@data$age[mymodel5@data$student == 4], yhat[mymodel5@data$student == 4], col = 4)

plot(mymodel5@data$age, yhat, type = "n")
lines(mymodel5@data$age[mymodel5@data$student == 10], yhat[mymodel5@data$student == 10], col = 1)
lines(mymodel5@data$age[mymodel5@data$student == 11], yhat[mymodel5@data$student == 11], col = 2)
lines(mymodel5@data$age[mymodel5@data$student == 12], yhat[mymodel5@data$student == 12], col = 3)
lines(mymodel5@data$age[mymodel5@data$student == 13], yhat[mymodel5@data$student == 13], col = 4)
lines(mymodel5@data$age[mymodel5@data$student == 14], yhat[mymodel5@data$student == 14], col = 4)

# Chapter learning outcomes . . . . . . . . . . . . . . . . . . . . . . .209

############################################################################
