## ---- include=FALSE------------------------------------------------------
library(knitr)
opts_chunk$set(out.extra='style="display:block; margin: auto"', fig.align="center")

## ---- eval = TRUE, warning = FALSE, fig.width = 7, fig.height = 5--------
################################################################################
### R code for
### "Global, Parameterwise and Joint Shrinkage Factor Estimation"
### written by Daniela Dunkler, March 2016
################################################################################

## works with R 3.2.3 & shrink-package 1.2.1

## load the packages
library("shrink")
library("survival")
library("mfp")
library("aod")                                                  # for wald-test

sessionInfo()

## ---- eval = TRUE, warning = FALSE, fig.width = 7, fig.height = 5--------
## Section 2.1: Deep Vein Thrombosis Study
## load data
data("deepvein")

# number of observations, events, median observation time, and names of
# variables
nrow(deepvein)

deepvein$status2 <- abs(deepvein$status-1)
survfit(Surv(time, status2) ~ 1, data = deepvein)

table(deepvein$status)
round(100*table(deepvein$status)/nrow(deepvein), 2)

sort(names(deepvein))

## ---- eval = TRUE, warning = FALSE, fig.width = 7, fig.height = 5--------
## Section 2.2: Breast Cancer Study
## load data
data("GBSG")


# number of observations, events, median observation time, and names of
# variables
nrow(GBSG)

table(GBSG$cens)
round(100*table(GBSG$cens)/nrow(GBSG), 2)

GBSG$cens2 <- abs(GBSG$cens-1)
# median observation time is given in days here
survfit(Surv(rfst, cens2) ~ 1, data = GBSG)
# median observation time in months
1645 / 30.5


sort(names(GBSG))

## ---- eval = TRUE, warning = FALSE, fig.width = 7, fig.height = 5-----------------
## Section 5.1: Deep Vein Thrombosis Study
# set the reference level for all categorical variables
deepvein$sex <- relevel(deepvein$sex, ref = "female")
deepvein$fiimut <- relevel(deepvein$fiimut, ref = "absent")
deepvein$fvleid <- relevel(deepvein$fvleid, ref = "absent")
deepvein$loc <- relevel(deepvein$loc, ref = "PE")
# contrasts(deepvein$loc)


## fit full model, and compute shrinkage factors (jackknife estimates and dfbeta
## approximations)
fitfull <- coxph(Surv(time, status) ~ log2ddim + bmi + durther + age +
                 sex + loc + fiimut + fvleid, data = deepvein, x = TRUE)
summary(fitfull)


## wald-test for categorical predictor "loc"
wald.test(Sigma = vcov(fitfull)[c("locdistal", "locproximal"),
          c("locdistal", "locproximal")],
          b = fitfull$coeff[c("locdistal", "locproximal")], Terms = 1:2)


shrink(fitfull, type = "global", method = "jackknife")



## apply backward elimination to full model*, and compute shrinkage factors to
## selected model
## * dummy variables of loc are jointly tested and selected
fitd <- step(fitfull, direction = "backward")
print(fitd)


## wald-test for categorical predictor "loc"
wald.test(Sigma = vcov(fitd)[c("locdistal", "locproximal"),
          c("locdistal", "locproximal")],
          b = fitd$coeff[c("locdistal", "locproximal")], Terms = 1:2)



shrink(fitd, type = "global", method = "jackknife")
pwsf <- shrink(fitd, type = "parameterwise", method = "jackknife")
pwsf

sqrt(diag(pwsf$ShrinkageFactorsVCOV))

shrink(fitd, type = "parameterwise", method = "jackknife",
       join = list(c("locdistal", "locproximal")))


## dummy variables of loc are separately tested and selected
## generate dummy variables, fit full model and then apply backward elimination
deepvein <- cbind(deepvein, model.matrix( ~ loc, data = deepvein)[, -1])
fitfull2 <- coxph(Surv(time, status) ~ log2ddim + bmi + durther + age +
                  sex + locdistal + locproximal + fiimut + fvleid,
                  data = deepvein, x = TRUE)  # fitfull2 is identical to fitfull

fitd2 <- step(fitfull2, direction = "backward")
summary(fitd2)

shrink(fitd2, type = "global", method = "jackknife")
shrink(fitd2, type = "parameterwise", method = "jackknife")



## Table 2
t2_jack <- c(shrink(fitd, type = "global", method = "jackknife")$ShrinkageFactors,
             pwsf$ShrinkageFactors,
             shrink(fitd, type = "parameterwise", method = "jackknife",
                    join = list(c("locdistal", "locproximal")))$ShrinkageFactors,
             system.time(shrink(fitd, type = "parameterwise", method = "jackknife"))[1])

t2_dfbeta <- c(shrink(fitd, type = "global", method = "dfbeta")$ShrinkageFactors,
               shrink(fitd, type = "parameterwise", method = "dfbeta")$ShrinkageFactors,
               shrink(fitd, type = "parameterwise", method = "dfbeta",
                      join = list(c("locdistal", "locproximal")))$ShrinkageFactors,
               system.time(shrink(fitd, type = "parameterwise", method = "dfbeta"))[1])

t2 <- cbind(Jackknife = t2_jack, DFBETA = t2_dfbeta)

Table2 <- cbind(t2, round((t2[, 2] - t2[, 1]) / t2[, 1] * 100, 1))
Table2[, 1:2] <- round(Table2[, 1:2], 4)
dimnames(Table2)[[2]][3] <- "Relative difference"
Table2 <- rbind(Table2[1,], rep(NA, 3), Table2[2:5,], rep(NA, 3), Table2[6:8,],
                rep(NA, 3), Table2[10,])
dimnames(Table2)[[1]][c(1:2, 7, 10, 12)] <- c("Global shrinkage",
          "Parameterwise shrinkage", "Joint shrinkage", "loc", "Computing time")
Table2

# max. difference
max(abs(Table2[1:10, 1] - Table2[1:10, 2]), na.rm=TRUE)



## Figure 1: Data simulated from deep vein thrombosis study: Comparison of Jackknife and
##           DFBETA-approximated global shrinkage factors for selected and unselected
##           models.

deepvein0 <- subset(deepvein, status == 0)
deepvein1 <- subset(deepvein, status == 1)

n0 <- nrow(deepvein0)
n1 <- nrow(deepvein1)

ratio <- 0.2                                     # based on n1 / n0
## Note: Running this code is time-consuming, since 2 * B * S data sets are analyzed.
## size <- seq(from = 200, to = 2000, length.out = 19)      # based on nrow(deepvein)
# B <- 100
## You may want to reduce B and length.out of size: e.g.
## size <- seq(from = 200, to = 2000, length.out = 4)
## B <- 3
#
# S <- length(size)
# 
# shrinkGJ <- shrinkGD <- shrinkGJsel <- shrinkGDsel <- matrix(NA, nrow = B, ncol = S)
# set.seed(954681)
# 
# for (s in 1:S) {
#   cat("\n", s)
#   for (b in 1:B) {
#     cat(".")
#     data <- rbind(deepvein0[sample(x = 1:n0, size = size[s]*(1-ratio), replace = TRUE),],
#                   deepvein1[sample(x = 1:n1, size = size[s]*ratio, replace = TRUE),])
#     
#     fitfull <- coxph(Surv(time, status) ~ log2ddim + bmi + durther + age + sex + loc +
#                        fiimut + fvleid, data = data, x = TRUE)
#     fitsel <- step(fitfull, direction = "backward", trace = 0)
#     
#     if (!is.null(fitsel$coefficients)) {                     # no null model selected
#       shrinkGJsel[b, s] <- shrink(fitsel, type = "global", method = "jackknife")$ShrinkageFactors
#       shrinkGDsel[b, s] <- shrink(fitsel, type = "global", method = "dfbeta")$ShrinkageFactors
#     }
#     
#     
#     fit <- coxph(Surv(time, status) ~ log2ddim + sex + loc, data = data, x = TRUE)
#     
#     shrinkGJ[b, s] <- shrink(fit, type = "global", method = "jackknife")$ShrinkageFactors
#     shrinkGD[b, s] <- shrink(fit, type = "global", method = "dfbeta")$ShrinkageFactors
#   }
# }
## In some smaller data sets (and especially when shrinkage factors are estimated with the
## Jackknife method) there might be issues with convergence in individual data sets.
## However, coxph will issue a warning.
# 
# 
# shrinkGDselm <- apply(shrinkGDsel, 2, median)
# shrinkGJselm <- apply(shrinkGJsel, 2, median)
# 
# cbind(n = size, Diff_J_D_sel = shrinkGJselm - shrinkGDselm)
# 
# shrinkGDm <- apply(shrinkGD, 2, median)
# shrinkGJm <- apply(shrinkGJ, 2, median)
# 
# cbind(n = size, Diff_J_D = shrinkGJm - shrinkGDm)
# 
# xrange <- c(0.5, 1)
## xrange <- range(shrinkGDselm, shrinkGJselm, shrinkGDm, shrinkGJm)

## pdf("compJvsD.pdf", width = 6, height = 4)
# par(oma = c(0, 0, 0.5, 0.5), mar = c(3.5, 4, 0, 0))
# plot(range(size), xrange, type = "n", las = 1, xlab = "", ylab = "", xaxt = "n")
# axis(1, at = size)
# mtext(side = 1, text = "Sample size", line = 2.3)
# mtext(side = 2, text = "Global shrinkage factor", line = 2.8)
# 
# points(size, shrinkGDselm, pch = 3, col = "darkolivegreen4", cex = 1.5)
# points(size, shrinkGJselm, pch = 1, col = "darkgoldenrod3", cex = 1.3)
# lines(size, shrinkGDselm, lty = 2, col = "darkolivegreen4", lwd = 2)
# lines(size, shrinkGJselm, lty = 3, col = "darkgoldenrod3", lwd = 2)
# 
# points(size, shrinkGDm, pch = 3, col = "firebrick3", cex = 1.5)
# points(size, shrinkGJm, pch = 1, col = "dodgerblue3", cex = 1.3)
# lines(size, shrinkGDm, lty = 2, col = "firebrick3", lwd = 2)
# lines(size, shrinkGJm, lty = 3, col = "dodgerblue3", lwd = 2)
# 
# legend(x = 1600, y = 0.62, legend = c("Jackknife", "DFBETA"), lwd = 2, title = "unselected",
#        col = c("dodgerblue3", "firebrick3"), inset = 0.05, bty = "n", pch = c(1, 3),
#        text.col = c("dodgerblue3", "firebrick3"), title.col = "black")
# 
# legend(x = 1000, y = 0.62, legend = c("Jackknife", "DFBETA"), lwd = 2, title = "selected",
#        col = c("darkgoldenrod3", "chartreuse4"), inset = 0.05, bty = "n", pch = c(1, 3),
#        text.col = c("darkgoldenrod3", "chartreuse4"), title.col = "black")
## dev.off()

## ---- eval = TRUE, warning = FALSE, fig.width = 7, fig.height = 5-----------------
## Section 5.2: Breast Cancer Study
## define predictors as suggested by Sauerbrei (1999, Applied Statistics)
GBSG$enodes <- exp(-0.12*GBSG$posnodal)


# create ordinal dummy-coded variable tumgrad for grades 1 to 3
contrasts(GBSG$tumgrad) <- matrix(c(0, 1, 1, 0, 0, 1), ncol = 2, byrow = FALSE,
                                  dimnames = list(1:3, c("tumgrad1", "tumgrad2")))
contrasts(GBSG$tumgrad)
# for nicer variable names use default dimnames
contrasts(GBSG$tumgrad) <- matrix(c(0, 1, 1, 0, 0, 1), ncol = 2, byrow = FALSE)


## # model selection (backward elimination) and estimation
fitg <- mfp(Surv(rfst, cens) ~ fp(age) + fp(prm) + fp(esm) + fp(tumsize) +
              fp(enodes) + tumgrad + menostat + strata(htreat),
            family = cox, data = GBSG, alpha = 0.05, select = 0.05)
print(fitg)


## Table 3
## (p-values unconditional on the selected degree and power for prm, age, and
## enodes = fitg$pvalues)
varorder <- c("age.1", "age.2", "prm.1", "enodes.1", "tumgrad1.1")

t3_beta_j <- c(NA, round(fitg$coefficients, 3)[varorder])
t3_df <- c(fitg$df.initial[7, ], NA, NA, fitg$df.initial[c(2, 1, 3), ])
t3_p <- c(round(fitg$pvalues[7, "p.null"], 3), NA, NA,
          round(fitg$pvalues[c(2, 1, 3), "p.null"], 3))
t3_pwsf0 <- shrink(fitg, type = "parameterwise", method = "jackknife")
t3_pwsf <- round(c(NA, t3_pwsf0$ShrinkageFactors[varorder]), 3)
t3_pwsfse <- round(c(NA, sqrt(diag(t3_pwsf0$ShrinkageFactorsVCOV))[varorder]), 3)
t3_jointsf0 <- shrink(fitg, type = "parameterwise", method = "jackknife",
                      join=list(c("age.1", "age.2")))
t3_jointsf <- round(c(t3_jointsf0$ShrinkageFactors[varorder[1]], NA, NA,
                      t3_jointsf0$ShrinkageFactors[varorder[-c(1:2)]]), 3)
t3_jointsfse <- round(c(sqrt(diag(t3_jointsf0$ShrinkageFactorsVCOV))["join.age.1"],
                        NA, NA,
                        sqrt(diag(t3_jointsf0$ShrinkageFactorsVCOV))[varorder[-c(1:2)]]), 3)

Table3 <- cbind("beta_j" = t3_beta_j, "df" = t3_df, "p value" = t3_p,
                "PWSF" = t3_pwsf, "PWSF SE" = t3_pwsfse,
                "Joint shrinkage" = t3_jointsf,
                "Joint shrinkage SE" = t3_jointsfse)
dimnames(Table3)[[1]][1] <- "age"
Table3



## compute shrinkage factors
globalsf <- shrink(fitg, type = "global", method = "jackknife")
globalsf
sqrt(globalsf$ShrinkageFactorsVCOV)

pwsf <- shrink(fitg, type = "parameterwise", method = "jackknife")
pwsf
round(sqrt(diag(pwsf$ShrinkageFactorsVCOV)), 3)
round(cov2cor(pwsf$ShrinkageFactorsVCOV), 3)

jointsf <- shrink(fitg, type = "parameterwise", method = "jackknife",
                  join=list(c("age.1", "age.2")))
jointsf
round(sqrt(diag(jointsf$ShrinkageFactorsVCOV)), 3)



## Figure 2: selected model: Log hazard relative to 50 years
# refit model fitg explicitly including the two dummy variables of
# tumgrad
GBSG <- cbind(GBSG, model.matrix(~tumgrad, data = GBSG)[, -1])

fitg <- mfp(Surv(rfst, cens) ~ fp(age) + fp(prm) + fp(esm) + fp(tumsize) +
              fp(enodes) + tumgrad1 + tumgrad2 + menostat + strata(htreat),
            family = cox, data = GBSG, alpha = 0.05, select = 0.05)

globalsf <- shrink(fitg, type = "global", method = "jackknife")
pwsf <- shrink(fitg, type = "parameterwise", method = "jackknife")
jointsf <- shrink(fitg, type = "parameterwise", method = "jackknife",
                  join=list(c("age.1", "age.2")))

# define new data for which prediction is requested
# newdatref is the new data set for the reference observation
age <- 30:75
newdat <- data.frame(age = age, enodes = 1, prm = 1, tumgrad1 = 0,
                     htreat = factor(1))
newdatref <- data.frame(age = 50, enodes = 1, prm = 1, tumgrad1 = 0,
                        htreat = factor(1))


xaxis <- seq(from = min(age), to = max(age), by = 5)


# pdf("plotgbsg.pdf", width = 6, height = 6)
par(fig = c(0, 1, 0, 0.3), new = FALSE, oma = c(0, 0, 0.5, 0.5),
    mar = c(3.5, 4, 0, 0))
hist(GBSG$age[(GBSG$age >= min(age)) & (GBSG$age <= max(age))],
     breaks = seq(from = min(age), to = max(age), by = 2.5), main = "",
     xlim = range(age), xaxs = "r", yaxs = "r", xlab = "", ylab = "",
     axes = FALSE, col = "gray88")
axis(side = 1, at = xaxis)
axis(side = 2, at = seq(from = 0, to = 60, by = 20), las = 2)
mtext(side = 1, text = "Age", line = 2.3)
mtext(side = 2, text = "Frequency", line = 2.8)
box()

par(fig = c(0,1,0.3,1), new = TRUE, oma = c(0, 0, 0.5, 0.5),
    mar = c(1, 4, 0, 0))
plot(range(age), c(-0.05, 0.8), xlab = "", ylab = "", type = "n", yaxs = "i",
     xaxt = "n", yaxt = "n")
axis(side = 1, at = xaxis, labels = rep("", times = length(xaxis)))
axis(side = 2, at = seq(from = 0, to = 0.7, by = 0.1), las = 2,
     labels = seq(from = 0, to = 0.7, by = 0.1))
mtext(side = 2, text = "Log hazard relative to 50 years", line = 2.8)
lines(x = age, y = predict(fitg, newdata = newdat, type = "lp") -
        predict(fitg, newdata = newdatref, type = "lp"), lty = 1, col = "gray25",
      lwd = 2)
lines(x = age, y = predict(globalsf, newdata = newdat, type = "lp") -
        predict(globalsf, newdata = newdatref, type = "lp"), lty = 5,
      col = "forestgreen", lwd = 2)
lines(x = age, y = predict(pwsf, newdata = newdat, type="lp") -
        predict(pwsf, newdata = newdatref, type = "lp"), lty = 2,
      col = "dodgerblue3", lwd = 2)
lines(x = age, y = predict(jointsf, newdata = newdat, type = "lp") -
        predict(jointsf, newdata = newdatref, type = "lp"), lty = 4,
      col = "firebrick3", lwd = 2)

legend("topright", inset = 0.02, bty = "n", lty = c(1, 5, 4, 2),
       title = "SHRINKAGE",
       legend = c("No", "Global", "Joint", "Parameterwise"),
       col = c("gray25", "forestgreen", "firebrick3", "dodgerblue3"), lwd = 2,
       text.col = c("gray25", "forestgreen", "firebrick3", "dodgerblue3"))
# dev.off()

