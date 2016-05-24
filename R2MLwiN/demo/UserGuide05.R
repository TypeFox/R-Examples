############################################################################
#     MLwiN User Manual
#
# 5   Graphical Procedures for Exploring the Model . . . . . . . . . . . .65
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


# 5.1 Displaying multiple graphs . . . . . . . . . . . . . . . . . . . . .65

data(tutorial, package = "R2MLwiN")

(mymodel1 <- runMLwiN(normexam ~ 1 + standlrt + (1 | school) + (1 | student), estoptions = list(resi.store = TRUE), 
  data = tutorial))

u0 <- mymodel1@residual$lev_2_resi_est_Intercept
u0se <- sqrt(mymodel1@residual$lev_2_resi_var_Intercept)

u0rank <- rank(u0)
u0rankhi <- u0 + u0se
u0ranklo <- u0 - u0se
u0rankno <- order(u0rank)

plot(1:65, u0[u0rankno], ylim = c(-1, 1), pch = 15, xlab = "Rank", ylab = "u0 residual estimate")
points(1:65, u0rankhi[u0rankno], pch = 24, bg = "grey")
points(1:65, u0ranklo[u0rankno], pch = 25, bg = "grey")
for (i in 1:65) {
  lines(rep(i, 2), c(u0ranklo[u0rankno[i]], u0rankhi[u0rankno[i]]))
}

plot(mymodel1@data$standlrt, mymodel1@data$normexam, asp = 1)

# 5.2 Highlighting in graphs . . . . . . . . . . . . . . . . . . . . . . .68

xb <- predict(mymodel1)

xbu <- xb + u0[mymodel1@data$school]

pred <- as.data.frame(cbind(mymodel1@data$school, mymodel1@data$standlrt, xb, xbu)[order(mymodel1@data$school, mymodel1@data$standlrt), 
  ])

colnames(pred) <- c("school", "standlrt", "xb", "xbu")

xyplot(xbu ~ standlrt, type = "l", group = school, data = pred)

xyplot(xbu ~ standlrt, panel = function(x, y, subscripts) {
  panel.xyplot(x, y, type = "l", groups = pred$school, subscripts = subscripts)
  panel.xyplot(pred$standlrt, pred$xb, type = "l", lwd = 3, color = "black")
}, data = pred)

c(unique(mymodel1@data$school)[u0rank == 65], u0rank[u0rank == 65])

sch53 <- which(levels(as.factor(mymodel1@data$school)) == 53)

plot(1:65, u0[u0rankno], ylim = c(-1, 1), pch = 15, xlab = "Rank", ylab = "u0 residual estimate")
points(1:65, u0rankhi[u0rankno], pch = 24, bg = "grey")
points(1:65, u0ranklo[u0rankno], pch = 25, bg = "grey")
for (i in 1:65) lines(rep(i, 2), c(u0ranklo[u0rankno[i]], u0rankhi[u0rankno[i]]))
points(x = which(u0rankno == sch53), y = u0[u0rankno[which(u0rankno == sch53)]], pch = 22, bg = 2)
legend(5, 1, "School 53", pch = 22, pt.bg = 2, col = 2)

plot(mymodel1@data$standlrt, xbu, type = "n")
for (i in 1:65) {
  lines(mymodel1@data$standlrt[mymodel1@data$school == i], xbu[mymodel1@data$school == i], col = "blue")
}
lines(mymodel1@data$standlrt, xb, col = 1, lwd = 3)
lines(mymodel1@data$standlrt[mymodel1@data$school == 53], xbu[mymodel1@data$school == 53], col = "red")
legend(-3, 2, "School 53", lty = 1, col = "red")

plot(mymodel1@data$standlrt, mymodel1@data$normexam, type = "p")
points(mymodel1@data$standlrt[mymodel1@data$school == 53], mymodel1@data$normexam[mymodel1@data$school == 53], col = "red")
legend(-3, 3, "School 53", lty = 1, col = "red")

schid <- as.vector(by(mymodel1@data$school, mymodel1@data$school, function(x) x[1]))
schsize <- as.vector(by(mymodel1@data$school, mymodel1@data$school, length))

cbind(schid, schsize, u0rank)[u0rank >= 28 & u0rank <= 32, ]

sch48 <- which(levels(as.factor(mymodel1@data$school)) == 48)

plot(1:65, u0[u0rankno], ylim = c(-1, 1), pch = 15, xlab = "Rank", ylab = "u0 residual estimate")
points(1:65, u0rankhi[u0rankno], pch = 24, bg = "grey")
points(1:65, u0ranklo[u0rankno], pch = 25, bg = "grey")
for (i in 1:65) {
  lines(rep(i, 2), c(u0ranklo[u0rankno[i]], u0rankhi[u0rankno[i]]))
}
points(x = which(u0rankno == sch53), y = u0[u0rankno[which(u0rankno == sch53)]], pch = 22, bg = 2)
points(x = which(u0rankno == sch48), y = u0[u0rankno[which(u0rankno == sch48)]], pch = 22, bg = 3)
legend(5, 1, c("School 53", "School 48"), pch = 22, pt.bg = c(2, 3), col = c(2, 3))

plot(mymodel1@data$standlrt, xbu, type = "n")
for (i in 1:65) {
  lines(mymodel1@data$standlrt[mymodel1@data$school == i], xbu[mymodel1@data$school == i], col = "blue")
}
lines(mymodel1@data$standlrt, xb, col = 1, lwd = 3)
lines(mymodel1@data$standlrt[mymodel1@data$school == 53], xbu[mymodel1@data$school == 53], col = "red")
lines(mymodel1@data$standlrt[mymodel1@data$school == 48], xbu[mymodel1@data$school == 48], col = "green")
legend(-3, 2, c("The average school", "School 53", "School 48"), lty = 1, col = c("black", "red", "green"))

plot(mymodel1@data$standlrt, mymodel1@data$normexam, type = "p")
points(mymodel1@data$standlrt[mymodel1@data$school == 53], mymodel1@data$normexam[mymodel1@data$school == 53], col = "red")
points(mymodel1@data$standlrt[mymodel1@data$school == 48], mymodel1@data$normexam[mymodel1@data$school == 48], col = "green")
legend(-3, 3, c("School 53", "School 48"), lty = 1, col = c("red", "green"))

cbind(schid, u0rank)[u0rank == 1, ]

sch59 <- which(levels(as.factor(mymodel1@data$school)) == 59)

plot(1:65, u0[u0rankno], ylim = c(-1, 1), pch = 15, xlab = "Rank", ylab = "u0 residual estimate")
points(1:65, u0rankhi[u0rankno], pch = 24, bg = "grey")
points(1:65, u0ranklo[u0rankno], pch = 25, bg = "grey")
for (i in 1:65) lines(rep(i, 2), c(u0ranklo[u0rankno[i]], u0rankhi[u0rankno[i]]))
points(x = which(u0rankno == sch53), y = u0[u0rankno[which(u0rankno == sch53)]], pch = 22, bg = 2)
points(x = which(u0rankno == sch59), y = u0[u0rankno[which(u0rankno == sch59)]], pch = 22, bg = 3)
legend(5, 1, c("School 53", "School 59"), pch = 22, pt.bg = c(2, 3), col = c(2, 3))

plot(mymodel1@data$standlrt, xbu, type = "n")
for (i in 1:65) {
  lines(mymodel1@data$standlrt[mymodel1@data$school == i], xbu[mymodel1@data$school == i], col = "blue")
}
lines(mymodel1@data$standlrt, xb, col = 1, lwd = 3)
lines(mymodel1@data$standlrt[mymodel1@data$school == 53], xbu[mymodel1@data$school == 53], col = "red")
lines(mymodel1@data$standlrt[mymodel1@data$school == 59], xbu[mymodel1@data$school == 59], col = "green")
legend(-3, 2, c("The average school", "School 53", "School 59"), lty = 1, col = c("black", "red", "green"))

plot(mymodel1@data$standlrt, mymodel1@data$normexam, type = "p")
points(mymodel1@data$standlrt[mymodel1@data$school == 53], mymodel1@data$normexam[mymodel1@data$school == 53], col = "red")
points(mymodel1@data$standlrt[mymodel1@data$school == 59], mymodel1@data$normexam[mymodel1@data$school == 59], col = "green")
legend(-3, 3, c("School 53", "School 59"), lty = 1, col = c("red", "green"))

(mymodel2 <- runMLwiN(normexam ~ 1 + standlrt + (1 + standlrt | school) + (1 | student), estoptions = list(resi.store = TRUE, 
  resioptions = c("estimates", "sampling")), data = tutorial))

xb <- predict(mymodel2)

u0 <- mymodel2@residual$lev_2_resi_est_Intercept
u1 <- mymodel2@residual$lev_2_resi_est_standlrt

xbu <- xb + u0[mymodel2@data$school] + u1[mymodel2@data$school] * mymodel2@data$standlrt

plot(u1, u0, xlab = "Slope", ylab = "Intercept")
points(u1[schid == 53], u0[schid == 53], col = "red")
points(u1[schid == 59], u0[schid == 59], col = "green")
legend(0.2, -0.5, c("School 53", "School 59"), pch = 22, pt.bg = c("red", "green"), col = c("red", "green"))

plot(mymodel2@data$standlrt, xbu, type = "n")
for (i in 1:65) {
  lines(mymodel2@data$standlrt[mymodel2@data$school == i], xbu[mymodel2@data$school == i], col = "blue")
}
lines(mymodel2@data$standlrt, xb, col = 1, lwd = 3)
lines(mymodel2@data$standlrt[mymodel2@data$school == 53], xbu[mymodel2@data$school == 53], col = "red")
lines(mymodel2@data$standlrt[mymodel2@data$school == 59], xbu[mymodel2@data$school == 59], col = "green")
legend(-3, 2, c("School 53", "School 59"), lty = 1, col = c("red", "green"))

plot(mymodel2@data$standlrt, mymodel2@data$normexam, type = "p")
points(mymodel2@data$standlrt[mymodel2@data$school == 53], mymodel2@data$normexam[mymodel2@data$school == 53], col = "red")
points(mymodel2@data$standlrt[mymodel2@data$school == 59], mymodel2@data$normexam[mymodel2@data$school == 59], col = "green")
legend(-3, 3, c("School 53", "School 59"), lty = 1, col = c("red", "green"))

u0var <- mymodel2@residual$lev_2_resi_var_Intercept
u0u1cov <- mymodel2@residual$lev_2_resi_cov_standlrt_Intercep
u1var <- mymodel2@residual$lev_2_resi_var_standlrt

xbu_lo <- xbu - 1.96 * sqrt(u0var[mymodel2@data$school] + 2 * u0u1cov[mymodel2@data$school] * mymodel2@data$standlrt + 
  u1var[mymodel2@data$school] * mymodel2@data$standlrt^2)
xbu_hi <- xbu + 1.96 * sqrt(u0var[mymodel2@data$school] + 2 * u0u1cov[mymodel2@data$school] * mymodel2@data$standlrt + 
  u1var[mymodel2@data$school] * mymodel2@data$standlrt^2)

plotdata <- as.data.frame(cbind(mymodel2@data$standlrt, mymodel2@data$school, xb, xbu, xbu_lo, xbu_hi)[order(mymodel2@data$standlrt), 
  ])
colnames(plotdata) <- c("standlrt", "school", "xb", "xbu", "xbu_lo", "xbu_hi")

plot(plotdata$standlrt, plotdata$xb, xlim = c(-4, 3), ylim = c(-2.5, 3), type = "l")
lines(plotdata$standlrt[plotdata$school == 53], plotdata$xbu[plotdata$school == 53], col = "red")
lines(plotdata$standlrt[plotdata$school == 53], plotdata$xbu_lo[plotdata$school == 53], lty = 3, col = "red")
lines(plotdata$standlrt[plotdata$school == 53], plotdata$xbu_hi[plotdata$school == 53], lty = 3, col = "red")
lines(plotdata$standlrt[plotdata$school == 59], plotdata$xbu[plotdata$school == 59], col = "green")
lines(plotdata$standlrt[plotdata$school == 59], plotdata$xbu_lo[plotdata$school == 59], lty = 3, col = "green")
lines(plotdata$standlrt[plotdata$school == 59], plotdata$xbu_hi[plotdata$school == 59], lty = 3, col = "green")
legend(-4, 3, c("School 53", "School 59"), lty = 1, col = c("red", "green"))

#     Chapter learning outcomes . . . . . . . . . . . . . . . . . . . . . 77

############################################################################
