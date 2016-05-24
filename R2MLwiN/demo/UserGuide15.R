############################################################################
#     MLwiN User Manual
#
# 15  Diagnostics for Multilevel Models . . . . . . . . . . . . . . . . .225
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


# 15.1 Introduction . . . . . . . . . . . . . . . . . . . . . . . . . . .225

data(diag1, package = "R2MLwiN")
summary(diag1)

(mymodel1 <- runMLwiN(n_ilea ~ 1 + n_vrq + (1 + n_vrq | school) + (1 | pupil), estoptions = list(resi.store = TRUE, 
  resioptions = c("standardised", "leverage", "influence", "deletion")), data = diag1))

xb <- predict(mymodel1)

u0 <- mymodel1@residual$lev_2_resi_est_Intercept
u1 <- mymodel1@residual$lev_2_resi_est_n_vrq

yhat <- xb + u0[mymodel1@data$school] + u1[mymodel1@data$school] * mymodel1@data$n_vrq

plot(mymodel1@data$n_vrq, yhat, type = "n")
for (i in 1:18) {
  lines(mymodel1@data$n_vrq[mymodel1@data$school == i], yhat[mymodel1@data$school == i], col = 1)
}
lines(mymodel1@data$n_vrq[mymodel1@data$school == 17], yhat[mymodel1@data$school == 17], col = 2, lwd = 3)

yhatold <- yhat

plot(mymodel1@data$n_vrq[mymodel1@data$school != 17], mymodel1@data$n_ilea[mymodel1@data$school != 17], type = "p")
points(mymodel1@data$n_vrq[mymodel1@data$school == 17], mymodel1@data$n_ilea[mymodel1@data$school == 17], col = 2, 
  cex = 2)

u0se <- mymodel1@residual$lev_2_resi_se_Intercept
u1se <- mymodel1@residual$lev_2_resi_se_n_vrq

u0rank <- rank(u0)
u0rankhi <- u0 + u0se
u0ranklo <- u0 - u0se
u0rankno <- order(u0rank)

u1rank <- rank(u1)
u1rankhi <- u1 + u1se
u1ranklo <- u1 - u1se
u1rankno <- order(u1rank)

sch17 <- which(levels(as.factor(mymodel1@data$school)) == 17)

plot(1:18, u0[u0rankno], ylim = c(-0.3, 0.3), pch = 15, xlab = "Rank", ylab = "u0 residual estimate")
points(1:18, u0rankhi[u0rankno], pch = 24, bg = "grey")
points(1:18, u0ranklo[u0rankno], pch = 25, bg = "grey")
for (i in 1:18) lines(rep(i, 2), c(u0ranklo[u0rankno[i]], u0rankhi[u0rankno[i]]))
points(x = which(u0rankno == sch17), y = u0[u0rankno[which(u0rankno == sch17)]], pch = 22, bg = 2, cex = 2)

plot(1:18, u1[u1rankno], ylim = c(-0.3, 0.3), pch = 15, xlab = "Rank", ylab = "u1 residual estimate")
points(1:18, u1rankhi[u1rankno], pch = 24, bg = "grey")
points(1:18, u1ranklo[u1rankno], pch = 25, bg = "grey")
for (i in 1:18) lines(rep(i, 2), c(u1ranklo[u1rankno[i]], u1rankhi[u1rankno[i]]))
points(x = which(u1rankno == sch17), y = u1[u1rankno[which(u1rankno == sch17)]], pch = 22, bg = 2, cex = 2)

plot(u0[-17], u1[-17], ylim = c(-0.3, 0.3), xlim = c(-0.1, 0.2), type = "p", xlab = "u0 residual estimate", ylab = "u1 residual estimate")
points(u0[17], u1[17], bg = 2, col = 2, cex = 2)

# 15.2 Diagnostics plotting: Deletion residuals, influence and leverage . 231

hist(u0, breaks = seq(min(u0) - 0.01, max(u0) + 0.01, 0.01))

u0std <- mymodel1@residual$lev_2_std_resi_est_Intercept
hist(u0std, breaks = seq(min(u0std) - 0.1, max(u0std) + 0.1, 0.1))

u0lev <- mymodel1@residual$lev_2_resi_leverage_Intercept
hist(u0lev, breaks = seq(min(u0lev) - 0.01, max(u0lev) + 0.01, 0.01))

u0inf <- mymodel1@residual$lev_2_resi_influence_Intercept
hist(u0inf, breaks = seq(min(u0inf) - 0.025, max(u0inf) + 0.025, 0.025))

u0del <- mymodel1@residual$lev_2_resi_deletion_Intercept
hist(u0del, breaks = seq(min(u0del) - 0.1, max(u0del) + 0.1, 0.1))

plot(u0std[-17], u0lev[-17], ylim = c(0.15, 0.4), xlim = c(-0.2, 2.2), type = "p", xlab = "u0 standardised residual", 
  ylab = "u0 leverage residual")
points(u0std[17], u0lev[17], bg = 2, col = 2, cex = 2)

hist(u1, breaks = seq(min(u1) - 0.01, max(u1) + 0.01, 0.01))

u1std <- mymodel1@residual$lev_2_std_resi_est_n_vrq
hist(u1std, breaks = seq(min(u1std) - 0.1, max(u1std) + 0.1, 0.1))

u1lev <- mymodel1@residual$lev_2_resi_leverage_n_vrq
hist(u1lev, breaks = seq(min(u1lev) - 0.01, max(u1lev) + 0.01, 0.01))

u1inf <- mymodel1@residual$lev_2_resi_influence_n_vrq
hist(u1inf, breaks = seq(min(u1inf) - 0.025, max(u1inf) + 0.025, 0.025))

u1del <- mymodel1@residual$lev_2_resi_deletion_n_vrq
hist(u1del, breaks = seq(min(u1del) - 0.1, max(u1del) + 0.1, 0.1))

plot(u1std[-17], u1lev[-17], ylim = c(0.1, 0.35), xlim = c(-1.4, 2.4), type = "p", xlab = "u1 standardised residual", 
  ylab = "u1 leverage residual")
points(u1std[17], u1lev[17], bg = 2, col = 2, cex = 2)

e0 <- mymodel1@residual$lev_1_resi_est_Intercept
e0rank <- rank(e0)
e0std <- (e0 - mean(e0))/sd(e0)
e0uniform <- e0rank/(length(e0rank) + 1)
e0nscore <- qnorm(e0uniform)

plot(e0nscore[mymodel1@data$school != 17], e0std[mymodel1@data$school != 17], ylim = c(-4, 5), xlim = c(-4, 4), type = "p", 
  xlab = "e0nscore", ylab = "e0std")
points(e0nscore[mymodel1@data$school == 17], e0std[mymodel1@data$school == 17], bg = 2, col = 2, cex = 2)

diag1$pupilnumber <- unlist(by(diag1$school, diag1$school, function(x) 1:length(x)))

diag1[order(e0)[1], c("school", "pupil", "pupilnumber")]

diag1$s17p22 <- as.integer(diag1$school == 17 & diag1$pupilnumber == 22)

(mymodel2 <- runMLwiN(n_ilea ~ 1 + n_vrq + s17p22 + (1 + n_vrq | school) + (1 | pupil), estoptions = list(resi.store = TRUE, 
  startval = list(FP.b = mymodel1@FP, FP.v = mymodel1@FP.cov, RP.b = mymodel1@RP, RP.v = mymodel1@RP.cov)), data = diag1))

u0 <- mymodel2@residual$lev_2_resi_est_Intercept
u1 <- mymodel2@residual$lev_2_resi_est_n_vrq

u0se <- sqrt(mymodel2@residual$lev_2_resi_var_Intercept)
u1se <- sqrt(mymodel2@residual$lev_2_resi_var_n_vrq)

u0rank <- rank(u0)
u0rankhi <- u0 + u0se
u0ranklo <- u0 - u0se
u0rankno <- order(u0rank)

u1rank <- rank(u1)
u1rankhi <- u1 + u1se
u1ranklo <- u1 - u1se
u1rankno <- order(u1rank)

plot(1:18, u0[u0rankno], ylim = c(-0.3, 0.3), pch = 15, xlab = "Rank", ylab = "u0 residual estimate")
points(1:18, u0rankhi[u0rankno], pch = 24, bg = "grey")
points(1:18, u0ranklo[u0rankno], pch = 25, bg = "grey")
for (i in 1:18) lines(rep(i, 2), c(u0ranklo[u0rankno[i]], u0rankhi[u0rankno[i]]))
points(x = which(u0rankno == sch17), y = u0[u0rankno[which(u0rankno == sch17)]], pch = 22, bg = 2, cex = 2)

plot(1:18, u1[u1rankno], ylim = c(-0.3, 0.3), pch = 15, xlab = "Rank", ylab = "u1 residual estimate")
points(1:18, u1rankhi[u1rankno], pch = 24, bg = "grey")
points(1:18, u1ranklo[u1rankno], pch = 25, bg = "grey")
for (i in 1:18) lines(rep(i, 2), c(u1ranklo[u1rankno[i]], u1rankhi[u1rankno[i]]))
points(x = which(u1rankno == sch17), y = u1[u1rankno[which(u1rankno == sch17)]], pch = 22, bg = 2, cex = 2)

diag1$s17 <- as.integer(diag1$school == 17)
diag1$s17Xn_vrq <- diag1$s17 * diag1$n_vrq

(mymodel3 <- runMLwiN(n_ilea ~ 1 + n_vrq + s17p22 + s17 + s17Xn_vrq + (1 + n_vrq | school) + (1 | pupil), estoptions = list(startval = list(FP.b = mymodel2@FP, 
  FP.v = mymodel2@FP.cov, RP.b = mymodel2@RP, RP.v = mymodel2@RP.cov)), data = diag1))

(mymodel4 <- runMLwiN(n_ilea ~ 1 + n_vrq + s17p22 + s17 + (1 + n_vrq | school) + (1 | pupil), estoptions = list(resi.store = TRUE, 
  startval = list(FP.b = mymodel3@FP, FP.v = mymodel3@FP.cov, RP.b = mymodel3@RP, RP.v = mymodel3@RP.cov)), data = diag1))

xb <- predict(mymodel4)

u0 <- mymodel4@residual$lev_2_resi_est_Intercept
u1 <- mymodel4@residual$lev_2_resi_est_n_vrq

yhat <- xb + u0[mymodel4@data$school] + u1[mymodel4@data$school] * mymodel4@data$n_vrq

plot(mymodel4@data$n_vrq, yhat, type = "n")
for (i in 1:18) {
  lines(mymodel4@data$n_vrq[mymodel4@data$school == i], yhatold[mymodel4@data$school == i], col = 1)
}
lines(mymodel4@data$n_vrq[mymodel4@data$school == 17], yhatold[mymodel4@data$school == 17], col = 2, lwd = 3)

plot(mymodel4@data$n_vrq, yhat, type = "n")
for (i in 1:18) {
  points(mymodel4@data$n_vrq[mymodel4@data$school == i], yhat[mymodel4@data$school == i], col = 1)
}
points(mymodel4@data$n_vrq[mymodel4@data$school == 17], yhat[mymodel4@data$school == 17], col = 2, lwd = 3)

# 15.3 A general approach to data exploration . . . . . . . . . . . . . .240

# Chapter learning outcomes . . . . . . . . . . . . . . . . . . . . . . .240

############################################################################
