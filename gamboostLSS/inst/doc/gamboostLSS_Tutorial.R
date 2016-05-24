### R code from vignette source 'gamboostLSS_Tutorial.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: init
###################################################
options(prompt = "R> ", continue = "+  ", width = 76, useFancyQuotes = FALSE)
repos <- "http://cran.at.r-project.org"
## check if a current version of gamboostLSS is installed
suppressWarnings(pd <- packageDescription("gamboostLSS"))
if (all(is.na(pd))){
    install.packages("gamboostLSS", repos = repos)
    pd <- packageDescription("gamboostLSS")
} else {
    if (compareVersion(pd$Version, "1.2-0") < 0){
        warning(sQuote("gamboostLSS"), " (1.2-0 or newer) is being installed!")
        install.packages("gamboostLSS", repos = repos)
    }
}

## for the exact reproduction of the plots R2BayesX >= 0.3.2 is needed
suppressWarnings(if (!require("R2BayesX"))
                     install.packages("R2BayesX"))

## remove (possibly) altered versions of "india" from working environment
suppressWarnings(rm("india"))
require("gamboostLSS")

## make graphics directory if not existing
if (!file.exists("graphics"))
    dir.create("graphics")
if (!file.exists("cvrisk"))
    dir.create("cvrisk")


###################################################
### code chunk number 2: load_package
###################################################
library("gamboostLSS")


###################################################
### code chunk number 3: simulate_data
###################################################
set.seed(1907)
n <- 150
x1 <- rnorm(n)
x2 <- rnorm(n)
x3 <- rnorm(n)
toydata <- data.frame(x1 = x1, x2 = x2, x3 = x3)
toydata$y <- rnorm(n, mean = 1 + 2 * x1 - x2,
                      sd = exp(0.5 - 0.25 * x1 + 0.5 * x3))


###################################################
### code chunk number 4: modelfitting
###################################################
lmLSS <- glmboostLSS(y ~ x1 + x2 + x3, data = toydata)


###################################################
### code chunk number 5: coefficients2
###################################################
coef(lmLSS, off2int = TRUE)


###################################################
### code chunk number 6: coef_paths (eval = FALSE)
###################################################
## par(mfrow = c(1, 2), mar = c(4, 4, 2, 5))
## plot(lmLSS, off2int = TRUE)


###################################################
### code chunk number 7: plot_coef_paths
###################################################
par(mfrow = c(1, 2), mar = c(4, 4, 2, 5))
plot(lmLSS, off2int = TRUE)


###################################################
### code chunk number 8: gamboostLSS_Tutorial.Rnw:308-310
###################################################
muFit <- fitted(lmLSS, parameter = "mu")
rbind(muFit, truth = 1 + 2 * x1 - x2)[, 1:5]


###################################################
### code chunk number 9: gamboostLSS_Tutorial.Rnw:315-317
###################################################
sigmaFit <- fitted(lmLSS, parameter = "sigma", type = "response")[, 1]
rbind(sigmaFit, truth = exp(0.5 - 0.25 * x1 + 0.5 * x3))[, 1:5]


###################################################
### code chunk number 10: load_package_and_data
###################################################
data("india")
data("india.bnd")
names(india)


###################################################
### code chunk number 11: india_stunting1
###################################################
library("R2BayesX")
## plot mean
FUN <- mean
india_agg <- data.frame(mcdist = names(tapply(india$stunting, india$mcdist, FUN,
                                              na.rm = TRUE)),
                        stunting = tapply(india$stunting, india$mcdist, FUN,
                                          na.rm = TRUE))
par(mar = c(1, 0, 2, 0))
plotmap(india.bnd, india_agg, pos = "bottomright",
        range = c(-5.1, 1.50), main = "Mean", mar.min = NULL)


###################################################
### code chunk number 12: india_stunting2
###################################################
## plot sd
FUN <- sd
india_agg <- data.frame(mcdist = names(tapply(india$stunting, india$mcdist, FUN,
                                              na.rm = TRUE)),
                        stunting = tapply(india$stunting, india$mcdist, FUN,
                                          na.rm = TRUE))
## remove missing values (= no variation)
india_agg <- india_agg[complete.cases(india_agg),]
par(mar = c(1, 0, 2, 0))
plotmap(india.bnd, india_agg, pos = "bottomright",
        range = c(0, 4.00), main = "Standard deviation")


###################################################
### code chunk number 13: summary
###################################################
out <- data.frame(Description = c("Stunting", "BMI (child)", "Age (child; months)",
                                 "BMI (mother)", "Age (mother; years)"),
                  Variable = paste0("\\code{", c(names(india)[1:5]),"}"),
                  "Min." = apply(india[, 1:5], 2, min),
                  "25\\% Quantile" = apply(india[, 1:5], 2, quantile, p = 0.25),
                  "Median" = apply(india[, 1:5], 2, median),
                  "Mean" = apply(india[, 1:5], 2, mean),
                  "75\\% Quantile" = apply(india[, 1:5], 2, quantile, p = 0.75),
                  "Max." = apply(india[, 1:5], 2, max))
names(out)[c(4, 7)] <- c("25\\% Qu.", "75\\% Qu.")
names(out)[1:2] <- ""
out[, 3:8] <- round(out[, 3:8], 2)
cat("\\toprule \n")
cat(paste(colnames(out), collapse = " & "), "\\\\ \n \\cmidrule{3-8}", "\n")
out <- apply(out, 1, function(x) cat(paste(x, collapse = " & "), " \\\\ \n"))
cat("\\bottomrule \n")


###################################################
### code chunk number 14: families
###################################################
str(GaussianLSS(), 1)


###################################################
### code chunk number 15: compute_neighborhood
###################################################
library("R2BayesX")
neighborhood <- bnd2gra(india.bnd)


###################################################
### code chunk number 16: options
###################################################
ctrl <- boost_control(trace = TRUE, mstop = c(mu = 1269, sigma = 84))


###################################################
### code chunk number 17: modeling (eval = FALSE)
###################################################
## mod_nonstab <- gamboostLSS(stunting ~ bbs(mage) + bbs(mbmi) +
##                              bbs(cage) + bbs(cbmi) +
##                              bmrf(mcdist, bnd = neighborhood),
##                            data = india,
##                            families = GaussianLSS(),
##                            control = ctrl)


###################################################
### code chunk number 18: run_modeling
###################################################
## abbreviate output to use less manuscript space
out <- capture.output(
mod_nonstab <- gamboostLSS(stunting ~ bbs(mage) + bbs(mbmi) +
                             bbs(cage) + bbs(cbmi) +
                             bmrf(mcdist, bnd = neighborhood),
                           data = india,
                           families = GaussianLSS(),
                           control = ctrl)
)
out <- out[c(1:5, 33:35)]
out[3:5] <- c("", "(...)", "")


###################################################
### code chunk number 19: print_modeling_results
###################################################
writeLines(out)


###################################################
### code chunk number 20: modeling2 (eval = FALSE)
###################################################
## mod <- gamboostLSS(stunting ~ bbs(mage) + bbs(mbmi) +
##                      bbs(cage) + bbs(cbmi) +
##                      bmrf(mcdist, bnd = neighborhood),
##                    data = india,
##                    families = GaussianLSS(stabilization = "MAD"),
##                    control = ctrl)


###################################################
### code chunk number 21: run_modeling2
###################################################
## abbreviate output to use less manuscript space
out <- capture.output(
mod <- gamboostLSS(stunting ~ bbs(mage) + bbs(mbmi) +
                     bbs(cage) + bbs(cbmi) +
                     bmrf(mcdist, bnd = neighborhood),
                   data = india,
                   families = GaussianLSS(stabilization = "MAD"),
                   control = ctrl)
)
out <- out[c(1:5, 33:35)]
out[3:5] <- c("", "(...)", "")


###################################################
### code chunk number 22: print_modeling2_results
###################################################
writeLines(out)


###################################################
### code chunk number 23: gamboostLSS_Tutorial.Rnw:908-910
###################################################
emp_risk <- risk(mod, merge = TRUE)
tail(emp_risk, n = 1)


###################################################
### code chunk number 24: gamboostLSS_Tutorial.Rnw:913-915
###################################################
emp_risk_nonstab <- risk(mod_nonstab, merge = TRUE)
tail(emp_risk_nonstab, n = 1)


###################################################
### code chunk number 25: gamboostLSS_Tutorial.Rnw:927-929
###################################################
weights <- sample(c(rep(1, 2000), rep(0, 2000)))
mod_subset <- update(mod, weights = weights, risk = "oobag")


###################################################
### code chunk number 26: gamboostLSS_Tutorial.Rnw:932-940 (eval = FALSE)
###################################################
## mod_subset <- gamboostLSS(stunting ~ bbs(mage) + bbs(mbmi) +
##                      bbs(cage) + bbs(cbmi) +
##                      bmrf(mcdist, bnd = neighborhood),
##                    data = india,
##                    weights = weights,
##                    families = GaussianLSS(),
##                    control = boost_control(mstop = c(mu = 1269, sigma = 84),
##                                            risk = "oobag"))


###################################################
### code chunk number 27: gamboostLSS_Tutorial.Rnw:943-945
###################################################
mod_nonstab_subset <- update(mod_nonstab,
                             weights = weights, risk = "oobag")


###################################################
### code chunk number 28: gamboostLSS_Tutorial.Rnw:949-951
###################################################
tail(risk(mod_subset, merge = TRUE), 1)
tail(risk(mod_nonstab_subset, merge = TRUE), 1)


###################################################
### code chunk number 29: setup_grid1
###################################################
grid <- make.grid(max = c(mu = 500, sigma = 500), min = 20,
                  length.out = 10, dense_mu_grid = FALSE)


###################################################
### code chunk number 30: setup_grid2_code (eval = FALSE)
###################################################
## densegrid <- make.grid(max = c(mu = 500, sigma = 500), min = 20,
##                        length.out = 10, dense_mu_grid = TRUE)
## plot(densegrid, pch = 20, cex = 0.2,
##      xlab = "Number of boosting iterations (mu)",
##      ylab = "Number of boosting iterations (sigma)")
## abline(0,1)
## points(grid, pch = 20, col = "red")


###################################################
### code chunk number 31: setup_grid2
###################################################
par(mar = c(4, 4, 0.1, 0.1))
densegrid <- make.grid(max = c(mu = 500, sigma = 500), min = 20,
                       length.out = 10, dense_mu_grid = TRUE)
plot(densegrid, pch = 20, cex = 0.2,
     xlab = "Number of boosting iterations (mu)",
     ylab = "Number of boosting iterations (sigma)")
abline(0,1)
points(grid, pch = 20, col = "red")


###################################################
### code chunk number 32: iterations
###################################################
par(mar = c(4, 4, 0.1, 0.1))
plot(1:30, type = "n", xlab = "Number of boosting iterations (mu)",
     ylab = "Number of boosting iterations (sigma)")
j <- 0
for (i in 0:29) {
    lines(c(i, i + 1), c(j, j), col = "red", lwd = 2)
    points(i, j, col = "red", pch = 20, cex = 1.3)
    if (j < 15) {
        lines(c(i + 1, i + 1), c(j, j + 1), col = "red", lwd = 2)
        points(i + 1, j, col = "red", pch = 20, cex = 1.3)
        j <- j + 1
    }
}
points(30, 15, col = "red", pch = 20, cex = 1.3)
lines(c(16, 16), c(15, 16), lwd = 2)
lines(c(16, 30), c(16, 16), lwd = 2)
points(c(16:30), rep(16, 15), col = "black", pch = 20, cex = 1.3)
abline(0,1)
lines(c(-1, 15), c(15, 15), lty = "dotted")
lines(c(15, 15), c(-1, 15), lty = "dotted")
legend("topleft", legend = c("mstop = c(mu = 30, sigma = 15)", "mstop = c(mu = 30, sigma = 16)"),
       lty = "solid", pch = 20, cex = 1.1, col = c("red", "black"), bty = "n")


###################################################
### code chunk number 33: cvrisk
###################################################
cores <- ifelse(grepl("linux|apple", R.Version()$platform), 2, 1)
if (!file.exists("cvrisk/cvr_india.Rda")) {
    set.seed(1907)
    folds <- cv(model.weights(mod), type = "subsampling")
    densegrid <- make.grid(max = c(mu = 5000, sigma = 500), min = 20,
                           length.out = 10, dense_mu_grid = TRUE)
    cvr <- cvrisk(mod, grid = densegrid, folds = folds, mc.cores = cores)
    save("cvr", file = "cvrisk/cvr_india.Rda", compress = "xz")
}


###################################################
### code chunk number 34: load_cvrisk
###################################################
load("cvrisk/cvr_india.Rda")


###################################################
### code chunk number 35: crossvalidation_code (eval = FALSE)
###################################################
## plot(cvr)


###################################################
### code chunk number 36: crossvalidation
###################################################
plot(cvr)


###################################################
### code chunk number 37: mstop
###################################################
mstop(cvr)


###################################################
### code chunk number 38: check_initial_assignement
###################################################
## for the purpose of the tutorial we started already with the optimal number of
## boosting steps. Check if this is really true:
if (!isTRUE(all.equal(mstop(mod), mstop(cvr))))
    warning("Check mstop(mod) throughout the manuscript.")


###################################################
### code chunk number 39: subset
###################################################
mstop(mod) <- mstop(cvr)


###################################################
### code chunk number 40: effects
###################################################
par(mfrow = c(2, 5))
plot(mod)


###################################################
### code chunk number 41: smooth_effects_code (eval = FALSE)
###################################################
## par(mfrow = c(2, 4), mar = c(5.1, 4.5, 4.1, 1.1))
## plot(mod, which = "bbs", type = "l")


###################################################
### code chunk number 42: smooth_effects
###################################################
par(cex.axis = 1.3, cex.lab = 1.3, mar = c(4, 5, 2, 1))
par(mfrow = c(2, 4), mar = c(5.1, 4.5, 4.1, 1.1))
plot(mod, which = "bbs", type = "l")


###################################################
### code chunk number 43: smooth_effects_2
###################################################
plot(mod, which = "bbs", parameter = "mu")


###################################################
### code chunk number 44: smooth_effects_3
###################################################
plot(mod, which = 1:4, parameter = 1)


###################################################
### code chunk number 45: marginal_prediction_plot_code
###################################################
plot(predint(mod, pi = c(0.8, 0.9), which = "cbmi"),
     lty = 1:3, lwd = 3, xlab =  "BMI (child)",
     ylab = "Stunting score")


###################################################
### code chunk number 46: greater_mumbai
###################################################
points(stunting ~ cbmi, data = india, pch = 20,
       col = rgb(1, 0, 0, 0.5), subset = mcdist == "381")


###################################################
### code chunk number 47: marginal_prediction_plot
###################################################
plot(predint(mod, pi = c(0.8, 0.9), which = "cbmi"),
     lty = 1:3, lwd = 3, xlab =  "BMI (child)",
     ylab = "Stunting score")
points(stunting ~ cbmi, data = india, pch = 20,
       col = rgb(1, 0, 0, 0.5), subset = mcdist == "381")


###################################################
### code chunk number 48: spatial_effects_code1 (eval = FALSE)
###################################################
## fitted_mu <- fitted(mod, parameter = "mu", which = "mcdist",
##                     type = "response")
## fitted_sigma <- fitted(mod, parameter = "sigma", which = "mcdist",
##                        type = "response")


###################################################
### code chunk number 49: spatial_effects_code2 (eval = FALSE)
###################################################
## fitted_mu <- tapply(fitted_mu, india$mcdist, FUN = mean)
## fitted_sigma <- tapply(fitted_sigma, india$mcdist, FUN = mean)
## plotdata <- data.frame(region = names(fitted_mu),
##                        mu = fitted_mu, sigma = fitted_sigma)
## par(mfrow = c(1, 2), mar = c(1, 0, 2, 0))
## plotmap(india.bnd, plotdata[, c(1, 2)], range = c(-0.62, 0.82),
##         main = "Mean", pos = "bottomright", mar.min = NULL)
## plotmap(india.bnd, plotdata[, c(1, 3)], range = c(0.75, 1.1),
##         main = "Standard deviation", pos = "bottomright", mar.min = NULL)


###################################################
### code chunk number 50: spatial_effects
###################################################
fitted_mu <- fitted(mod, parameter = "mu", which = "mcdist",
                    type = "response")
fitted_sigma <- fitted(mod, parameter = "sigma", which = "mcdist",
                       type = "response")
fitted_mu <- tapply(fitted_mu, india$mcdist, FUN = mean)
fitted_sigma <- tapply(fitted_sigma, india$mcdist, FUN = mean)
plotdata <- data.frame(region = names(fitted_mu),
                       mu = fitted_mu, sigma = fitted_sigma)
par(mfrow = c(1, 2), mar = c(1, 0, 2, 0))
plotmap(india.bnd, plotdata[, c(1, 2)], range = c(-0.62, 0.82),
        main = "Mean", pos = "bottomright", mar.min = NULL)
plotmap(india.bnd, plotdata[, c(1, 3)], range = c(0.75, 1.1),
        main = "Standard deviation", pos = "bottomright", mar.min = NULL)


