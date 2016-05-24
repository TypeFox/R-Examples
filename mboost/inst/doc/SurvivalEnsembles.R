### R code from vignette source 'SurvivalEnsembles.Rnw'

###################################################
### code chunk number 1: setup
###################################################
source("setup.R")
if (!require("TH.data"))
    stop("cannot attach package ", sQuote("TH.data"))
if (!require("rpart"))
    stop("cannot attach package ", sQuote("rpart"))
if (!require("survival"))
    stop("cannot attach package ", sQuote("survival"))
if (!require("party"))
    stop("cannot attach package ", sQuote("party"))

set.seed(290875)
CEX <- 0.85

### mean difference plots
mdplot <- function(obs, pred, main = "", ...) {

    m <- (obs + pred)/2
    d <- obs - pred
    plot(m, d, xlab = "(Observed + Predicted) / 2",
         ylab = "Observed - Predicted", main =
         main, cex.axis = CEX, cex.main = CEX, cex.lab = CEX, ...)
    abline(h = 0, lty = 3)
}


###################################################
### code chunk number 2: loaddata
###################################################
### load data. See `mboost/inst/readAML_Bullinger.R' for
### how the data were generated from the raw data.
load(file.path(path.package(package = "TH.data"), "rda", "AML_Bullinger.rda"))


###################################################
### code chunk number 3: AML-dpp
###################################################
### compute IPC weights
AMLw <- IPCweights(Surv(clinical$time, clinical$event))

### risk score
risk <- rep(0, nrow(clinical))
rlev <- levels(clinical[, "Cytogenetic.group"])
risk[clinical[, "Cytogenetic.group"] %in% rlev[c(7,8,4)]] <- "low"
risk[clinical[, "Cytogenetic.group"] %in% rlev[c(5, 9)]] <- "intermediate"
risk[clinical[, "Cytogenetic.group"] %in% rlev[-c(4,5, 7,8,9)]] <- "high"
risk <- as.factor(risk)

### set-up learning sample
AMLlearn <- cbind(clinical[, c("time", "Sex", "Age", "LDH", "WBC",
                        "FLT3.aberration.", "MLL.PTD", "Tx.Group.")],
               risk = risk,
               iexpressions[, colnames(iexpressions) %in% selgenes[["Clone.ID"]]])
cc <- complete.cases(AMLlearn)
AMLlearn <- AMLlearn[AMLw > 0 & cc,]
AMLw <- AMLw[AMLw > 0 & cc]


###################################################
### code chunk number 4: AML-RF
###################################################
### controls for tree growing
ctrl <- cforest_control(mincriterion = 0.1, mtry = 5, minsplit = 5, ntree = 250)

### fit random forest for censored data (warnings are OK here)
AMLrf <- cforest(I(log(time)) ~ ., data = AMLlearn, control = ctrl,
                 weights = AMLw)


###################################################
### code chunk number 5: AML-boost
###################################################
AMLl2b <- glmboost(I(log(time)) ~ ., data = AMLlearn, weights = AMLw,
                    control = boost_control(mstop = 5000))


###################################################
### code chunk number 6: AML-AIC
###################################################
### AIC criterion
plot(aic <- AIC(AMLl2b))


###################################################
### code chunk number 7: AML-fitted
###################################################
### restrict number of boosting iterations and inspect selected variables
AMLl2b <- AMLl2b[mstop(aic)]
cAML <- coef(AMLl2b)
cAML[abs(cAML) > 0]

### fitted values
AMLprf <- predict(AMLrf, newdata = AMLlearn)
AMLpb <- predict(AMLl2b, newdata = AMLlearn)


###################################################
### code chunk number 8: Figure1
###################################################
Mmod <- sum(AMLw * log(AMLlearn$time))/sum(AMLw )
par(mai = par("mai") * c(0.7, 0.8, 0.7, 0.6))
layout(matrix(1:4, ncol = 2))

mdplot(log(AMLlearn$time), AMLprf, main = "Random Forest",
       cex = AMLw / 4, ylim = c(-4, 4), xlim = c(0, 7))

plot(log(AMLlearn$time), AMLprf, cex = AMLw / 4,
       ylim = range(log(AMLlearn$time)), ylab = "Predicted", xlab = "Observed",
       main = "Random Forest",  cex.axis = CEX, cex.main = CEX, cex.lab = CEX)
abline(h = Mmod, lty = 2)

mdplot(log(AMLlearn$time), AMLpb,
        cex = AMLw / 4, main = "Boosting", ylim = c(-4, 4), xlim = c(0, 7))

plot(log(AMLlearn$time), AMLpb, cex = AMLw / 4,
       ylim = range(log(AMLlearn$time)), ylab = "Predicted", xlab = "Observed",
       main = "Boosting",  cex.axis = CEX, cex.main = CEX, cex.lab = CEX)
abline(h = Mmod, lty = 2)


###################################################
### code chunk number 9: GBSG2-dpp
###################################################
### attach data
data("GBSG2", package = "TH.data")

### IPC weights
GBSG2w <- IPCweights(Surv(GBSG2$time, GBSG2$cens))

### set-up learning sample
GBSG2learn <- cbind(GBSG2[,-which(names(GBSG2) %in% c("time", "cens"))],
               ltime = log(GBSG2$time))
n <- nrow(GBSG2learn)


###################################################
### code chunk number 10: GBSG2-models
###################################################
### linear model
LMmod <- lm(ltime ~ . , data = GBSG2learn, weights = GBSG2w)
LMerisk <- sum((GBSG2learn$ltime - predict(LMmod))^2*GBSG2w) / n

### regression tree
pos <- GBSG2w > 0
TRmod <- rpart(ltime ~ . , data = GBSG2learn, weights = GBSG2w, 
               subset = pos)
TRerisk <- sum((GBSG2learn$ltime[pos] - predict(TRmod))^2*GBSG2w[pos]) / n

### tree controls
ctrl <- cforest_control(mincriterion = qnorm(0.95), mtry = 5,
                      minsplit = 5, ntree = 100)

### fit random forest for censored data (warnings are OK here)
RFmod <- cforest(ltime ~ . , data = GBSG2learn, weights = GBSG2w,
                 control = ctrl)

### fit L2 boosting for censored data
L2Bmod <- glmboost(ltime ~ ., data = GBSG2learn, weights = GBSG2w,
                   control = boost_control(mstop = 250))

### with Huber loss function
L2BHubermod <- glmboost(ltime ~ ., data = GBSG2learn, weights = GBSG2w,
                        family = Huber(d = log(2)))


###################################################
### code chunk number 11: GBSG2-AIC
###################################################
plot(aic <- AIC(L2Bmod))


###################################################
### code chunk number 12: GBSG2-fitted
###################################################
GBSG2Hp <- predict(L2BHubermod, newdata = GBSG2learn)
L2Berisk <- sum((GBSG2learn$ltime - predict(L2Bmod, newdata = GBSG2learn))^2*GBSG2w) / n
RFerisk <- sum((GBSG2learn$ltime - predict(RFmod, newdata = GBSG2learn))^2*GBSG2w) / n


###################################################
### code chunk number 13: Figure3
###################################################
lim <- c(4,9)
mylwd <- 0.5
par(mai = par("mai") * c(0.7, 0.8, 0.7, 0.6))

layout(matrix(1:4, ncol = 2))
Mmod <- sum(GBSG2w * GBSG2learn$ltime)/sum(GBSG2w)

mdplot(GBSG2learn$ltime, predict(LMmod), cex = GBSG2w / 4, main = "Linear Model",
       ylim = c(-3, 3), xlim = c(5, 8))
mdplot(GBSG2learn$ltime, predict(TRmod), cex = GBSG2w / 4, main = "Tree",
       ylim = c(-3, 3), xlim = c(5, 8))
mdplot(GBSG2learn$ltime, predict(RFmod, newdata = GBSG2learn), cex = GBSG2w / 4,
       main = "Random Forest", ylim = c(-3, 3), xlim = c(5, 8))
mdplot(GBSG2learn$ltime, predict(L2Bmod, newdata = GBSG2learn), cex = GBSG2w / 4,
       main = "Boosting", ylim = c(-3, 3), xlim = c(5, 8))


###################################################
### code chunk number 14: Figure5
###################################################
RFpr <- predict(RFmod, newdata = GBSG2learn)
L2Bpr <- predict(L2Bmod, newdata = GBSG2learn)
ylim <- range(c(RFpr[GBSG2w > 0], L2Bpr[GBSG2w > 0]))
mydf <- 4
par(mai = par("mai") * c(0.7, 0.8, 0.4, 0.6))
layout(matrix(1:4, ncol = 2))
plot(GBSG2learn$pnodes, RFpr, cex = GBSG2w/4, xlim = c(0,40), lwd = mylwd,
     xlab = "Nr. positive lymph nodes", ylim = ylim, ylab =
     expression(hat(Y)), cex.axis = CEX, cex.main = CEX, cex.lab = CEX,)
lines(smooth.spline(GBSG2learn$pnodes, RFpr, GBSG2w/4, df = mydf))
plot(GBSG2learn$age, RFpr, cex = GBSG2w/4, xlab = "Age",
     ylab = expression(hat(Y)), ylim = ylim, lwd = mylwd, cex.axis = CEX,
     cex.main = CEX, cex.lab = CEX)
lines(smooth.spline(GBSG2learn$age, RFpr, GBSG2w/4, df = mydf))
plot(GBSG2learn$estrec, RFpr, cex = GBSG2w/4, xlab = "Estrogen receptor",
     ylab = expression(hat(Y)), ylim = ylim, lwd = mylwd, cex.axis = CEX,
     cex.main = CEX, cex.lab = CEX)
lines(smooth.spline(GBSG2learn$estrec, RFpr, GBSG2w/4, df = mydf))
indx <- which(GBSG2learn$progrec < 100)
plot(GBSG2learn$progrec[indx], RFpr[indx], cex = GBSG2w[indx]/4,
     xlab = "Progesterone receptor (< 100 fmol / l)",
     ylab = expression(hat(Y)), ylim = ylim, lwd = mylwd, cex.axis = CEX,
     cex.main = CEX, cex.lab = CEX)
lines(smooth.spline(GBSG2learn$progrec[indx], RFpr[indx], GBSG2w[indx]/4, df = mydf))


###################################################
### code chunk number 15: Figure6
###################################################
par(mai = par("mai") * c(0.7, 0.8, 0.4, 0.6))
layout(matrix(1:4, ncol = 2))
plot(GBSG2learn$pnodes, L2Bpr, cex = GBSG2w/4, xlim = c(0,40),
     ylab = expression(hat(Y)),
     xlab = "Nr. positive lymph nodes", ylim = ylim, lwd = mylwd, cex.axis =
     CEX, cex.main = CEX, cex.lab = CEX)
lines(smooth.spline(GBSG2learn$pnodes, L2Bpr, GBSG2w/4, df = mydf))
plot(GBSG2learn$age, L2Bpr, cex = GBSG2w/4, xlab = "Age",
     ylab = expression(hat(Y)),
     ylim = ylim, lwd = mylwd, cex.axis = CEX, cex.main = CEX, cex.lab =
     CEX)
lines(smooth.spline(GBSG2learn$age, L2Bpr, GBSG2w/4, df = mydf))
plot(GBSG2learn$estrec, L2Bpr, cex = GBSG2w/4, xlab = "Estrogen receptor",
     ylab = expression(hat(Y)), ylim = ylim, lwd = mylwd, cex.axis = CEX,
     cex.main = CEX, cex.lab = CEX)
lines(smooth.spline(GBSG2learn$estrec, L2Bpr, GBSG2w/4, df = mydf))
indx <- which(GBSG2learn$progrec < 100)
plot(GBSG2learn$progrec[indx], L2Bpr[indx], cex = GBSG2w[indx]/4,
     xlab = "Progesterone receptor (< 100 fmol / l)",
     ylab = expression(hat(Y)), ylim = ylim, lwd = mylwd, cex.axis = CEX,
     cex.main = CEX, cex.lab = CEX)
lines(smooth.spline(GBSG2learn$progrec[indx], L2Bpr[indx], GBSG2w[indx]/4, df = mydf))


###################################################
### code chunk number 16: Figure7
###################################################
Mmod <- sum(GBSG2w * GBSG2learn$ltime)/sum(GBSG2w)
par(mai = par("mai") * c(0.7, 0.8, 0.7, 0.6))
layout(matrix(1:4, ncol = 2))

yl <- range(c(GBSG2Hp[GBSG2w > 0], L2Bpr[GBSG2w > 0]))
mdplot(GBSG2learn$ltime, GBSG2Hp, main = "Huber Loss",
       cex = GBSG2w / 4, ylim = c(-3, 3), xlim = c(5, 8))

plot(GBSG2learn$ltime, GBSG2Hp, cex = GBSG2w / 4,
       xlim = range(GBSG2learn$ltime[GBSG2w > 0]), ylim = yl, ylab = "Predicted", xlab = "Observed",
       main = "Huber Loss",  cex.axis = CEX, cex.main = CEX, cex.lab = CEX)

mdplot(GBSG2learn$ltime, L2Bpr,
        cex = GBSG2w / 4, main = "Quadratic Loss", ylim = c(-3, 3), xlim = c(5, 8))

plot(GBSG2learn$ltime, L2Bpr, cex = GBSG2w / 4,
       xlim = range(GBSG2learn$ltime[GBSG2w > 0]), ylim = yl, ylab = "Predicted", xlab = "Observed",
       main = "Quadratic Loss",  cex.axis = CEX, cex.main = CEX, cex.lab = CEX)


