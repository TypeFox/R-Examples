### R code from vignette source 'mboost_illustrations.Rnw'

###################################################
### code chunk number 1: pkg-attach
###################################################
source("setup.R")


###################################################
### code chunk number 2: bodyfat-lm-fit
###################################################
bf_lm <- lm(DEXfat ~ hipcirc + kneebreadth + anthro3a, data = bodyfat)
coef(bf_lm)


###################################################
### code chunk number 3: bodyfat-glmboost-fit
###################################################
bf_glm <- glmboost(DEXfat ~ ., data = bodyfat, center = TRUE)


###################################################
### code chunk number 4: bodyfat-glmboost-coef
###################################################
coef(bf_glm)


###################################################
### code chunk number 5: bodyfat-oob-plot
###################################################
load(system.file("cache/bodyfat_benchmarks.rda", package = "mboost"))
aic <- AIC(bf_glm)
pdf("figures/bodyfat_glmboost-bodyfat-oob-plot.pdf", version = "1.4", width = 6, height = 10)
par(mai = par("mai") * c(1, 1, 0.5, 1))
mopt <- grid[which.min(colMeans(boob))]
layout(matrix(1:2, nrow = 2))
perfplot(boob, grid, ylab = "Out-of-bootstrap squared error",
    xlab = "Number of boosting iterations", alpha = 0.05)
abline(h = mean(boobrest), lty = 2)
lines(c(which.min(colMeans(boob)), which.min(colMeans(boob))),
      c(0, min(colMeans(boob))), lty = 2)
points(which.min(colMeans(boob)), min(colMeans(boob)))
plot(aic, ylim = c(3, 5.5))
dev.off()


###################################################
### code chunk number 6: bodyfat-glmboost-AIC
###################################################
mstop(aic <- AIC(bf_glm))


###################################################
### code chunk number 7: bodyfat-glmboost-coef
###################################################
coef(bf_glm[mstop(aic)])


###################################################
### code chunk number 8: bodyfat-glmboost-coef-count
###################################################
cf <- coef(bf_glm[mopt])
nsel <- length(cf)


###################################################
### code chunk number 9: bodyfat-pkg-attach
###################################################
source("setup.R")


###################################################
### code chunk number 10: bodyfat-gamboost-fit
###################################################
bf_gam <- gamboost(DEXfat ~ ., data = bodyfat, baselearner = "bss")


###################################################
### code chunk number 11: bodyfat-gamboost-prune
###################################################
mstop(aic <- AIC(bf_gam))


###################################################
### code chunk number 12: bodyfat-gamboost-plot
###################################################
bf_gam <- bf_gam[mstop(aic)]
layout(matrix(1:4, ncol = 2))
plot(bf_gam, which = c("hipcirc", "kneebreadth", "waistcirc", "anthro3b"))


###################################################
### code chunk number 13: bodyfat-pkg-attach
###################################################
source("setup.R")
library("splines")
indep <- names(bodyfat)[names(bodyfat) != "DEXfat"]
bsfm <- as.formula(paste("DEXfat ~ ", paste("bs(", indep, ")", collapse = " + ", sep = ""), sep = ""))


###################################################
### code chunk number 14: bodyfat-bs
###################################################
bsfm


###################################################
### code chunk number 15: bodyfat-fpboost-fit
###################################################
ctrl <- boost_control(mstop = 5000)
bf_bs <- glmboost(bsfm, data = bodyfat, control = ctrl)
mstop(aic <- AIC(bf_bs))


###################################################
### code chunk number 16: bodyfat-fpboost-plot
###################################################
layout(matrix(1:4, ncol = 2, byrow = TRUE))
par(mai = par("mai") * c(1, 1, 0.5, 1))
cf <- coef(bf_bs[mstop(aic)])
varorder <- c("hipcirc", "waistcirc", "kneebreadth", "anthro3b")
fpartial <- sapply(varorder, function(v) {
    rowSums(predict(bf_bs, which = v))
})
out <- sapply(varorder, function(i) {
    x <- bodyfat[,i]
    plot(sort(x), fpartial[order(x),i],  main = "", type = "b",
         xlab = i, ylab = expression(f[partial]), ylim = max(abs(fpartial)) * c(-1, 1))
    abline(h = 0, lty = 2, lwd = 0.5)
    })


###################################################
### code chunk number 17: pkg-attach
###################################################
source("setup.R")

### OK, this once required Biobase and is a very dirty hack...
data("Westbc", package = "TH.data")
westbc <- Westbc
exprs <- function(x) x$assay
pData <- function(x) x$pheno

p <- nrow(exprs(westbc))
n <- ncol(exprs(westbc))
ytab <- table(pData(westbc)$nodal.y)


###################################################
### code chunk number 18: west-data-attach
###################################################
### extract matrix of expression levels and binary response
x <- t(exprs(westbc))
y <- pData(westbc)$nodal.y


###################################################
### code chunk number 19: west-glmboost-ytrans
###################################################
### numeric 0/1 response variable
yfit <- as.numeric(y) - 1


###################################################
### code chunk number 20: west-glmboost-family
###################################################
### L2 boosting for classification with response in 0/1
### and binomial log-likelihood as loss function
### ATTENTION: use offset = 1/2 instead of 0!!!
rho <- function(y, f, w = 1) {
    p <- pmax(pmin(1 - 1e-5, f), 1e-5)
    -y * log(p) - (1 - y) * log(1 - p)
}
ngradient <- function(y, f, w = 1) y - f
offset <- function(y, w) weighted.mean(y, w)
L2fm <- Family(ngradient = ngradient,
               loss = rho, offset = offset)


###################################################
### code chunk number 21: west-glmboost-fit (eval = FALSE)
###################################################
## ### fit a linear model with initial mstop = 200 boosting iterations
## ctrl <- boost_control(mstop = 200)
## west_glm <- glmboost(x, yfit, family = L2fm, center = TRUE,
##                      control = ctrl)


###################################################
### code chunk number 22: west-glmboost-fit-lars
###################################################
### fit a linear model with initial mstop = 200 boosting iterations
### and record time
time <- system.time(west_glm <- glmboost(x, yfit, family = L2fm, center = TRUE,
                    control = boost_control(mstop = 200)))[1]
x <- t(exprs(westbc) - rowMeans(exprs(westbc)))


###################################################
### code chunk number 23: west-glmboost-AIC
###################################################
### evaluate AIC based on binomial log-likelihood for _all_ boosting
### iterations m = 1, ..., mstop = 200
aic <- AIC(west_glm, method = "classical")
### where should one stop? mstop = 108 or 107
mstop(aic)


###################################################
### code chunk number 24: west-glmboost-cf
###################################################
cf <- coef(west_glm[mstop(aic)], which = 1:ncol(x))
ps <- length(cf[abs(cf) > 0])


###################################################
### code chunk number 25: west-glmboost-plot
###################################################
### compute standard deviations of expression levels for each gene
sdx <- sqrt(rowSums((t(x) - colMeans(x))^2)/(nrow(x) - 1))
layout(matrix(1:2, ncol = 2))
cf <- cf * sdx
plot(sort(cf[abs(cf) > 0]), ylab = "Standardized coefficients")
abline(h = 0, lty = 2)
plot(aic, ylim = c(25, 40))


###################################################
### code chunk number 26: pkg-attach
###################################################
source("setup.R")
n <- sum(complete.cases(wpbc))
p <- ncol(wpbc) - 2


###################################################
### code chunk number 27: wpbc-glm-fit
###################################################
### remove missing values and time variable
cc <- complete.cases(wpbc)
wpbc2 <- wpbc[cc, colnames(wpbc) != "time"]
### fit logistic regression model
wpbc_step <- step(glm(status ~ ., data = wpbc2, family = binomial()), trace = 0)


###################################################
### code chunk number 28: wpbc-glm-aic
###################################################
logLik(wpbc_step)
AIC(wpbc_step)


###################################################
### code chunk number 29: wpbc-glmboost-fit
###################################################
### fit logistic regression model via gradient boosting
ctrl <- boost_control(mstop = 500)
wpbc_glm <- glmboost(status ~ ., data = wpbc2, family = Binomial(),
                     center = TRUE, control = ctrl)


###################################################
### code chunk number 30: wpbc-glmboost-AIC
###################################################
aic <- AIC(wpbc_glm, "classical")
aic


###################################################
### code chunk number 31: wpbc-glmboost-fit2
###################################################
### fit with new mstop
wpbc_glm <- wpbc_glm[mstop(aic)]
coef(wpbc_glm)[abs(coef(wpbc_glm)) > 0]


###################################################
### code chunk number 32: wpbc-gamboost-fit
###################################################
wpbc_gam <- gamboost(status ~ ., data = wpbc2, family = Binomial(), baselearner = "bss")
mopt <- mstop(aic <- AIC(wpbc_gam, "classical"))
aic


###################################################
### code chunk number 33: wpbc-gamboost-plot
###################################################
wpbc_gam <- wpbc_gam[mopt]
fpartial <- predict(wpbc_gam, which = 1:length(variable.names(wpbc_gam)))
layout(matrix(1:4, nrow = 2, byrow = TRUE))
par(mai = par("mai") * c(1, 1, 0.5, 1))
plot(wpbc_gam, which = rev(order(colMeans(abs(fpartial))))[1:4])


###################################################
### code chunk number 34: pkg-attach
###################################################
source("setup.R")


###################################################
### code chunk number 35: wpbc-glmboost-PIC
###################################################
library("survival")
### calculate IPC weights
censored <- wpbc$status == "R"
iw <- IPCweights(Surv(wpbc$time, censored))
wpbc3 <- wpbc[,names(wpbc) != "status"]


###################################################
### code chunk number 36: wpbc-glmboost-censored-fit
###################################################
ctrl <- boost_control(mstop = 500)
wpbc_surv <- glmboost(log(time) ~ ., data = wpbc3,
                  weights = iw, center = TRUE, control = ctrl)
mstop(aic <- AIC(wpbc_surv))
wpbc_surv <- wpbc_surv[mstop(aic)]


###################################################
### code chunk number 37: wpbc-glmboost-coef
###################################################
names(coef(wpbc_surv)[abs(coef(wpbc_surv)) > 0])


###################################################
### code chunk number 38: wpbc-glmboost-censored-fit
###################################################
plot(log(wpbc3$time), predict(wpbc_surv, newdata = wpbc3),
     cex = iw, ylim = c(0, 5), xlim = c(0, 5),
     xlab = "Time to recurrence (log-scale)",
     ylab = "Predicted time to recurrence")
abline(a = 0, b = 1, lty = 2, lwd = 0.5)


