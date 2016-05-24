### R code from vignette source 'fitting_models.rnw'

###################################################
### code chunk number 1: Setup
###################################################
options(width=65)                    
##if(!require(systemfit, quietly=TRUE)) install.packages("systemfit")
require(lattice, quietly=TRUE)
lattice.options(default.theme = canonical.theme(color = FALSE))


###################################################
### code chunk number 2: sampletrees
###################################################
trees <- data.frame(plot=factor(c(1, 1, 1, 2, 2, 2)),
                    dbh.cm=c(30, 32, 35, 30, 33, 35),
                    ht.m=c(25, 30, 40, 30, 40, 50))


###################################################
### code chunk number 3: fig-cs-1
###################################################
plot(trees$dbh.cm, trees$ht.m, pch=c(1, 19)[trees$plot], 
     xlab="Diameter (cm)", ylab="Height (m)")
abline(lm(ht.m ~ dbh.cm, data=trees), col="darkgrey")


###################################################
### code chunk number 4: fig-cs-2
###################################################
case.model.1 <- lm(ht.m ~ dbh.cm, data=trees)
plot(fitted(case.model.1), residuals(case.model.1), 
     ylab = "Residuals", xlab = "Fitted Values", 
     pch = c(1, 19)[trees$plot])
abline(h = 0, col = "darkgrey")


###################################################
### code chunk number 5: fig-cs-3
###################################################
case.model.2 <- lm(ht.m ~ dbh.cm*plot, data=trees)
plot(fitted(case.model.2), residuals(case.model.2), 
     ylab = "Residuals", xlab = "Fitted Values", 
     pch = c(1, 19)[trees$plot])
abline(h = 0, col = "darkgrey")


###################################################
### code chunk number 6: case-data
###################################################
par(mar=c(4,4,2,2))
plot(trees$dbh.cm, trees$ht.m, pch=c(1, 19)[trees$plot], 
     xlab="Diameter (cm)", ylab="Height (m)")
abline(lm(ht.m ~ dbh.cm, data=trees), col="darkgrey")


###################################################
### code chunk number 7: case-res-2
###################################################
par(mar=c(4,4,2,2))
case.model.1 <- lm(ht.m ~ dbh.cm, data=trees)
plot(fitted(case.model.1), residuals(case.model.1), 
     ylab = "Residuals", xlab = "Fitted Values", 
     pch = c(1, 19)[trees$plot])
abline(h = 0, col = "darkgrey")


###################################################
### code chunk number 8: case-res-3
###################################################
par(mar=c(4,4,2,2))
case.model.2 <- lm(ht.m ~ dbh.cm*plot, data=trees)
plot(fitted(case.model.2), residuals(case.model.2), 
     ylab = "Residuals", xlab = "Fitted Values", 
     pch = c(1, 19)[trees$plot])
abline(h = 0, col = "darkgrey")


###################################################
### code chunk number 9: example
###################################################
example <- data.frame(y = c(4.2, 4.8, 5.8, 1.2, 10.1, 
                            14.9, 15.9, 13.1),
                      x = c(1, 2, 3, 4, 1, 2, 3, 4),
                      group = factor(c(1, 1, 1, 1, 
                                       2, 2, 2, 2)))


###################################################
### code chunk number 10: fitting_models.rnw:517-518
###################################################
library(nlme)


###################################################
### code chunk number 11: fig-simple1
###################################################
par(las=1, mar=c(4,4,1,1))
colours <- c("red", "blue")
plot(y ~ x, data = example, 
     col = colours[group], 
     pch = as.numeric(group))


###################################################
### code chunk number 12: simple1
###################################################
par(las=1, mar=c(4,4,1,1))
colours <- c("red", "blue")
plot(y ~ x, data = example, 
     col = colours[group], 
     pch = as.numeric(group))


###################################################
### code chunk number 13: fitting_models.rnw:573-574
###################################################
basic.1 <- lm(y ~ x, data=example)


###################################################
### code chunk number 14: fitting_models.rnw:578-579
###################################################
coef(summary(basic.1))


###################################################
### code chunk number 15: fitting_models.rnw:600-601
###################################################
basic.2 <- lm(y ~ x + group, data=example)


###################################################
### code chunk number 16: fitting_models.rnw:621-622
###################################################
basic.3 <- lm(y ~ x * group, data=example)


###################################################
### code chunk number 17: fitting_models.rnw:635-636
###################################################
example.mixed <- groupedData(y ~ x | group, data=example)


###################################################
### code chunk number 18: basic.4
###################################################
basic.4 <- lme(y ~ x, 
               random = ~1 | group, 
               data = example.mixed)


###################################################
### code chunk number 19: fig-simple4
###################################################
plot(augPred(basic.4))


###################################################
### code chunk number 20: simple4
###################################################
print(
plot(augPred(basic.4))
)


###################################################
### code chunk number 21: capture4
###################################################
summary(basic.4)


###################################################
### code chunk number 22: capture4
###################################################
summ.4 <- capture.output(summary(basic.4))


###################################################
### code chunk number 23: capture4a
###################################################
cat(summ.4[1:4], sep="\n")


###################################################
### code chunk number 24: capture4b
###################################################
cat(summ.4[6:9], sep="\n")


###################################################
### code chunk number 25: capture4c
###################################################
cat(summ.4[11:17], sep="\n")


###################################################
### code chunk number 26: capture4e
###################################################
cat(summ.4[19:21], sep="\n")


###################################################
### code chunk number 27: capture4f
###################################################
cat(summ.4[23:24], sep="\n")


###################################################
### code chunk number 28: basic.5
###################################################
basic.5 <- lme(y ~ x, random = ~1 | group,
           weights = varIdent(form = ~1 | group),
           data = example.mixed)


###################################################
### code chunk number 29: capture5
###################################################
summary(basic.5)


###################################################
### code chunk number 30: capture5
###################################################
summ.5 <- capture.output(summary(basic.5))


###################################################
### code chunk number 31: capture5f
###################################################
cat(summ.5[11:16], sep="\n")


###################################################
### code chunk number 32: basic.6
###################################################
basic.6 <- lme(y ~ x, random = ~1 | group,
               weights = varIdent(form = ~1 | group),
               correlation = corAR1(),
               data = example.mixed)


###################################################
### code chunk number 33: capture6
###################################################
summary(basic.6)


###################################################
### code chunk number 34: capture6
###################################################
summ.6 <- capture.output(summary(basic.6))


###################################################
### code chunk number 35: capture6f
###################################################
cat(summ.6[11:15], sep="\n")


###################################################
### code chunk number 36: get.stage
###################################################
Stangle("../../ch2/sweave/stage.rnw")
source("stage.R")


###################################################
### code chunk number 37: fitting_models.rnw:868-870
###################################################
names(stage)
dim(stage)


###################################################
### code chunk number 38: fitting_models.rnw:882-883
###################################################
stage.old <- stage[stage$Decade == 0, ]


###################################################
### code chunk number 39: fitting_models.rnw:922-924
###################################################
stage.old <- groupedData(height.m ~ dbhib.cm | Forest.ID,
                         data = stage.old)


###################################################
### code chunk number 40: fitting_models.rnw:1018-1020
###################################################
stage.old <- groupedData(height.m ~ dbhib.cm | Forest.ID, 
                         data = stage[stage$Decade == 0,])


###################################################
### code chunk number 41: fitting_models.rnw:1026-1029
###################################################
hd.lme.1 <- lme(height.m ~ dbhib.cm, 
                random = ~1 | Forest.ID,
                data = stage.old)


###################################################
### code chunk number 42: fig-hd-lme-1a
###################################################
scatter.smooth(fitted(hd.lme.1, level=0), 
               stage.old$height.m,
               xlab = "Fitted Values (height, m.)",
               ylab = "Observed Values (height, m.)",  
               main = "Model Structure (I)")
abline(0, 1, col = "blue")


###################################################
### code chunk number 43: fig-hd-lme-1b
###################################################
scatter.smooth(fitted(hd.lme.1), 
               residuals(hd.lme.1, type="pearson"),
               main = "Model Structure (II)",
               xlab = "Fitted Values", 
               ylab = "Innermost Residuals")
abline(h = 0, col = "red")


###################################################
### code chunk number 44: fig-hd-lme-1c
###################################################
ref.forest <- ranef(hd.lme.1)[[1]]
ref.var.forest <- 
  tapply(residuals(hd.lme.1, type="pearson", level=1),
         stage.old$Forest.ID,  var)
qqnorm(ref.forest, main="Q-Q Norm: Forest Random Effects")
qqline(ref.forest, col="red")


###################################################
### code chunk number 45: fig-hd-lme-1d
###################################################
qqnorm(residuals(hd.lme.1, type="pearson"), 
       main="Q-Q Normal - Residuals")
qqline(residuals(hd.lme.1, type="pearson"), col="red")


###################################################
### code chunk number 46: fig-hd-lme-1e
###################################################
boxplot(residuals(hd.lme.1, type = "pearson", level = 1) ~ 
        stage.old$Forest.ID,
        ylab = "Innermost Residuals", 
        xlab = "National Forest",
        notch = TRUE, 
        varwidth = TRUE, 
        at = rank(ref.forest))
axis(3, labels = format(ref.forest, dig=2), 
     cex.axis = 0.8, at = rank(ref.forest))
abline(h = 0, col = "darkgreen")


###################################################
### code chunk number 47: fig-hd-lme-1f
###################################################
plot(ref.forest, ref.var.forest, 
     xlab = "Forest Random Effect",
     ylab = "Variance of within-Forest Residuals")
abline(lm(ref.var.forest ~ ref.forest), col="purple")


###################################################
### code chunk number 48: diag-lme-1
###################################################
opar <- par(mfrow = c(3, 2), mar = c(4, 4, 3, 1), las = 1, cex.axis = 0.9)
scatter.smooth(fitted(hd.lme.1, level=0), 
               stage.old$height.m,
               xlab = "Fitted Values (height, m.)",
               ylab = "Observed Values (height, m.)",  
               main = "Model Structure (I)")
abline(0, 1, col = "blue")
scatter.smooth(fitted(hd.lme.1), 
               residuals(hd.lme.1, type="pearson"),
               main = "Model Structure (II)",
               xlab = "Fitted Values", 
               ylab = "Innermost Residuals")
abline(h = 0, col = "red")
ref.forest <- ranef(hd.lme.1)[[1]]
ref.var.forest <- 
  tapply(residuals(hd.lme.1, type="pearson", level=1),
         stage.old$Forest.ID,  var)
qqnorm(ref.forest, main="Q-Q Norm: Forest Random Effects")
qqline(ref.forest, col="red")
qqnorm(residuals(hd.lme.1, type="pearson"), 
       main="Q-Q Normal - Residuals")
qqline(residuals(hd.lme.1, type="pearson"), col="red")
boxplot(residuals(hd.lme.1, type = "pearson", level = 1) ~ 
        stage.old$Forest.ID,
        ylab = "Innermost Residuals", 
        xlab = "National Forest",
        notch = TRUE, 
        varwidth = TRUE, 
        at = rank(ref.forest))
axis(3, labels = format(ref.forest, dig=2), 
     cex.axis = 0.8, at = rank(ref.forest))
abline(h = 0, col = "darkgreen")
plot(ref.forest, ref.var.forest, 
     xlab = "Forest Random Effect",
     ylab = "Variance of within-Forest Residuals")
abline(lm(ref.var.forest ~ ref.forest), col="purple")
par(opar)


###################################################
### code chunk number 49: capturelme1
###################################################
summary.1 <- capture.output(summary(hd.lme.1))


###################################################
### code chunk number 50: capturelme1 (eval = FALSE)
###################################################
## summary(hd.lme.1)


###################################################
### code chunk number 51: capturelme1a
###################################################
cat(summary.1[1:4], sep="\n")


###################################################
### code chunk number 52: capturelme1b
###################################################
cat(summary.1[6:9], sep="\n")


###################################################
### code chunk number 53: capturelme1c
###################################################
cat(summary.1[11:14], sep="\n")


###################################################
### code chunk number 54: capturelme1d
###################################################
cat(summary.1[15:17], sep="\n")


###################################################
### code chunk number 55: capturelme1e
###################################################
cat(summary.1[19:21], sep="\n")


###################################################
### code chunk number 56: capturelme1f
###################################################
cat(summary.1[23:24], sep="\n")


###################################################
### code chunk number 57: fitting_models.rnw:1336-1337
###################################################
anova(hd.lme.1)


###################################################
### code chunk number 58: fitting_models.rnw:1366-1368
###################################################
stage <- groupedData(height.m ~ dbhib.cm | Forest.ID/Tree.ID,
                     data = stage)


###################################################
### code chunk number 59: fitting_models.rnw:1424-1427
###################################################
hd.lme.3 <- lme(height.m ~ dbhib.cm,
                random = ~1 | Forest.ID/Tree.ID,
                data = stage)


###################################################
### code chunk number 60: hd-lme-3a
###################################################
opar <- par(mfrow = c(1, 3), mar = c(4, 4, 3, 1), las = 1,
            cex.axis = 0.9)
plot(fitted(hd.lme.3, level=0), stage$height.m,
     xlab = "Fitted Values", ylab = "Observed Values",
     main = "Model Structure (I)")
abline(0, 1, col = "gray")
scatter.smooth(fitted(hd.lme.3), 
               residuals(hd.lme.3, type="pearson"),
               main = "Model Structure (II)",
               xlab = "Fitted Values", 
               ylab = "Innermost Residuals")
abline(h = 0, col = "gray")
acf.resid <- ACF(hd.lme.3, resType = "normal")
plot(acf.resid$lag[acf.resid$lag < 10.5],
     acf.resid$ACF[acf.resid$lag < 10.5],
     type="b", main="Autocorrelation",
     xlab="Lag", ylab="Correlation")
stdv <- qnorm(1 - 0.01/2)/sqrt(attr(acf.resid, "n.used"))
lines(acf.resid$lag[acf.resid$lag < 10.5],
      stdv[acf.resid$lag < 10.5],
      col="darkgray")
lines(acf.resid$lag[acf.resid$lag < 10.5],
      -stdv[acf.resid$lag < 10.5],
      col="darkgray")
abline(0,0,col="gray")
par(opar)


###################################################
### code chunk number 61: hd-lme-3a-d
###################################################
opar <- par(mfrow = c(1, 3), mar = c(4, 4, 3, 1), las = 1,
            cex.axis = 0.9)
plot(fitted(hd.lme.3, level=0), stage$height.m,
     xlab = "Fitted Values", ylab = "Observed Values",
     main = "Model Structure (I)")
abline(0, 1, col = "gray")
scatter.smooth(fitted(hd.lme.3), 
               residuals(hd.lme.3, type="pearson"),
               main = "Model Structure (II)",
               xlab = "Fitted Values", 
               ylab = "Innermost Residuals")
abline(h = 0, col = "gray")
acf.resid <- ACF(hd.lme.3, resType = "normal")
plot(acf.resid$lag[acf.resid$lag < 10.5],
     acf.resid$ACF[acf.resid$lag < 10.5],
     type="b", main="Autocorrelation",
     xlab="Lag", ylab="Correlation")
stdv <- qnorm(1 - 0.01/2)/sqrt(attr(acf.resid, "n.used"))
lines(acf.resid$lag[acf.resid$lag < 10.5],
      stdv[acf.resid$lag < 10.5],
      col="darkgray")
lines(acf.resid$lag[acf.resid$lag < 10.5],
      -stdv[acf.resid$lag < 10.5],
      col="darkgray")
abline(0,0,col="gray")
par(opar)


###################################################
### code chunk number 62: hd-lme-3b
###################################################
opar <- par(mfrow = c(1, 3), mar = c(4, 4, 3, 1), las = 1,
        cex.axis = 0.9)
ref.forest <- ranef(hd.lme.3, level=1, standard=T)[[1]]
ref.tree <- ranef(hd.lme.3, level=2, standard=T)[[1]]
ref.tree.frame <- ranef(hd.lme.3, level=2, 
                        augFrame=TRUE, standard=TRUE)
ref.var.tree <- tapply(residuals(hd.lme.3, type="pearson", 
                                 level=1),
                       stage$Tree.ID,  var)
ref.var.forest <- tapply(ref.tree, 
                         ref.tree.frame$Forest, var)
qqnorm(ref.forest, main = "QQ plot: Forest")
qqline(ref.forest)
qqnorm(ref.tree, main = "QQ plot: Tree")
qqline(ref.tree)
qqnorm(residuals(hd.lme.3, type="pearson"), 
       main="QQ plot: Residuals")
qqline(residuals(hd.lme.3, type="pearson"), col="red")
par(opar)


###################################################
### code chunk number 63: hd-lme-3b-d
###################################################
opar <- par(mfrow = c(1, 3), mar = c(4, 4, 3, 1), las = 1,
        cex.axis = 0.9)
ref.forest <- ranef(hd.lme.3, level=1, standard=T)[[1]]
ref.tree <- ranef(hd.lme.3, level=2, standard=T)[[1]]
ref.tree.frame <- ranef(hd.lme.3, level=2, 
                        augFrame=TRUE, standard=TRUE)
ref.var.tree <- tapply(residuals(hd.lme.3, type="pearson", 
                                 level=1),
                       stage$Tree.ID,  var)
ref.var.forest <- tapply(ref.tree, 
                         ref.tree.frame$Forest, var)
qqnorm(ref.forest, main = "QQ plot: Forest")
qqline(ref.forest)
qqnorm(ref.tree, main = "QQ plot: Tree")
qqline(ref.tree)
qqnorm(residuals(hd.lme.3, type="pearson"), 
       main="QQ plot: Residuals")
qqline(residuals(hd.lme.3, type="pearson"), col="red")
par(opar)


###################################################
### code chunk number 64: hd-lme-3c
###################################################
opar <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), 
            las = 1, cex.axis = 0.9)
boxplot(ref.tree ~ ref.tree.frame$Forest,
        ylab = "Tree Effects", xlab = "National Forest",
        notch= TRUE, varwidth = TRUE, at = rank(ref.forest))
axis(3, labels=format(ref.forest, dig=2), cex.axis=0.8,
     at=rank(ref.forest))
abline(h=0, col="darkgreen")
boxplot(residuals(hd.lme.3, type="pearson", level = 1) ~ 
        stage$Tree.ID,
        ylab = "Innermost Residuals", xlab = "Tree",
        notch = TRUE, varwidth = TRUE, at=rank(ref.tree))
axis(3, labels=format(ref.tree, dig=2), cex.axis=0.8,
     at=rank(ref.tree))
abline(h=0, col="darkgreen")
plot(ref.forest, ref.var.forest, xlab="Forest Random Effect",
     ylab="Variance of within-Forest Residuals")
abline(lm(ref.var.forest ~ ref.forest))
plot(ref.tree, ref.var.tree, xlab="Tree Random Effect",
     ylab="Variance of within-Tree Residuals")
abline(lm(ref.var.forest ~ ref.forest))

par(opar)


###################################################
### code chunk number 65: hd-lme-3c-d
###################################################
opar <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), 
            las = 1, cex.axis = 0.9)
boxplot(ref.tree ~ ref.tree.frame$Forest,
        ylab = "Tree Effects", xlab = "National Forest",
        notch= TRUE, varwidth = TRUE, at = rank(ref.forest))
axis(3, labels=format(ref.forest, dig=2), cex.axis=0.8,
     at=rank(ref.forest))
abline(h=0, col="darkgreen")
boxplot(residuals(hd.lme.3, type="pearson", level = 1) ~ 
        stage$Tree.ID,
        ylab = "Innermost Residuals", xlab = "Tree",
        notch = TRUE, varwidth = TRUE, at=rank(ref.tree))
axis(3, labels=format(ref.tree, dig=2), cex.axis=0.8,
     at=rank(ref.tree))
abline(h=0, col="darkgreen")
plot(ref.forest, ref.var.forest, xlab="Forest Random Effect",
     ylab="Variance of within-Forest Residuals")
abline(lm(ref.var.forest ~ ref.forest))
plot(ref.tree, ref.var.tree, xlab="Tree Random Effect",
     ylab="Variance of within-Tree Residuals")
abline(lm(ref.var.forest ~ ref.forest))

par(opar)


###################################################
### code chunk number 66: fitting_models.rnw:1628-1636
###################################################
library(lattice)
trees.in.forests <- 
  with(stage, aggregate(x = list(measures = height.m),
                        by = list(tree = Tree.ID, 
                                  forest = Forest.ID), 
                        FUN = length))
panel.order <- 
  rank(as.numeric(as.character(trees.in.forests$tree)))


###################################################
### code chunk number 67: hd.tree
###################################################
plot(augPred(hd.lme.3), 
     index.cond = list(panel.order),
     strip = strip.custom(par.strip.text = list(cex = 0.5)))


###################################################
### code chunk number 68: hd-tree
###################################################
print(
plot(augPred(hd.lme.3), 
     index.cond = list(panel.order),
     strip = strip.custom(par.strip.text = list(cex = 0.5)))
      )


###################################################
### code chunk number 69: hd.lme.4
###################################################
hd.lme.4 <- lme(height.m ~ dbhib.cm + I(dbhib.cm^2), 
                random = list( ~ 1 | Forest.ID, 
                               ~ dbhib.cm | Tree.ID),
                control = lmeControl(maxIter = 500, 
                                     msMaxIter = 500), 
                data = stage)


###################################################
### code chunk number 70: hd-lme-4
###################################################
opar <- par(mfrow = c(1, 3), mar = c(4, 4, 3, 1), las = 1,
        cex.axis = 0.9)
  plot(fitted(hd.lme.4, level=0), stage$height.m,
      xlab = "Fitted Values", ylab = "Observed Values",
      main = "Model Structure (I)")
  abline(0, 1, col = "gray")
  scatter.smooth(fitted(hd.lme.4), 
                 residuals(hd.lme.4, type="pearson"),
                 main = "Model Structure (II)",
                 xlab = "Fitted Values", 
                 ylab = "Innermost Residuals")
  abline(0, 0, col = "gray")
acf.resid <- ACF(hd.lme.4, resType = "n")
plot(acf.resid$lag[acf.resid$lag < 10.5],
     acf.resid$ACF[acf.resid$lag < 10.5],
     type="b", main="Autocorrelation",
     xlab="Lag", ylab="Correlation")
stdv <- qnorm(1 - 0.01/2)/sqrt(attr(acf.resid, "n.used"))
lines(acf.resid$lag[acf.resid$lag < 10.5],
      stdv[acf.resid$lag < 10.5],
      col="darkgray")
lines(acf.resid$lag[acf.resid$lag < 10.5],
      -stdv[acf.resid$lag < 10.5],
      col="darkgray")
abline(0,0,col="gray")
par(opar)


###################################################
### code chunk number 71: hd-lme-4-d
###################################################
opar <- par(mfrow = c(1, 3), mar = c(4, 4, 3, 1), las = 1,
        cex.axis = 0.9)
  plot(fitted(hd.lme.4, level=0), stage$height.m,
      xlab = "Fitted Values", ylab = "Observed Values",
      main = "Model Structure (I)")
  abline(0, 1, col = "gray")
  scatter.smooth(fitted(hd.lme.4), 
                 residuals(hd.lme.4, type="pearson"),
                 main = "Model Structure (II)",
                 xlab = "Fitted Values", 
                 ylab = "Innermost Residuals")
  abline(0, 0, col = "gray")
acf.resid <- ACF(hd.lme.4, resType = "n")
plot(acf.resid$lag[acf.resid$lag < 10.5],
     acf.resid$ACF[acf.resid$lag < 10.5],
     type="b", main="Autocorrelation",
     xlab="Lag", ylab="Correlation")
stdv <- qnorm(1 - 0.01/2)/sqrt(attr(acf.resid, "n.used"))
lines(acf.resid$lag[acf.resid$lag < 10.5],
      stdv[acf.resid$lag < 10.5],
      col="darkgray")
lines(acf.resid$lag[acf.resid$lag < 10.5],
      -stdv[acf.resid$lag < 10.5],
      col="darkgray")
abline(0,0,col="gray")
par(opar)


###################################################
### code chunk number 72: hd-lme-5
###################################################
hd.lme.5 <- update(hd.lme.4, correlation = corCAR1())


###################################################
### code chunk number 73: hd-lme-5-d
###################################################
opar <- par(mfrow = c(1, 3), mar = c(4, 4, 3, 1), las = 1,
        cex.axis = 0.9)
  plot(fitted(hd.lme.5, level=0), stage$height.m,
      xlab = "Fitted Values", ylab = "Observed Values",
      main = "Model Structure (I)")
  abline(0, 1, col = "gray")
  scatter.smooth(fitted(hd.lme.5), residuals(hd.lme.5, type="pearson"),
      main = "Model Structure (II)",
      xlab = "Fitted Values", ylab = "Innermost Residuals")
  abline(0, 0, col = "gray")
acf.resid <- ACF(hd.lme.5, resType = "n")
plot(acf.resid$lag[acf.resid$lag < 10.5],
     acf.resid$ACF[acf.resid$lag < 10.5],
     type="b", main="Autocorrelation",
     xlab="Lag", ylab="Correlation")
stdv <- qnorm(1 - 0.01/2)/sqrt(attr(acf.resid, "n.used"))
lines(acf.resid$lag[acf.resid$lag < 10.5],
      stdv[acf.resid$lag < 10.5],
      col="darkgray")
lines(acf.resid$lag[acf.resid$lag < 10.5],
      -stdv[acf.resid$lag < 10.5],
      col="darkgray")
abline(0,0,col="gray")
par(opar)


###################################################
### code chunk number 74: hab5a
###################################################
plot(hd.lme.5, resid(.) ~ fitted(.) | Hab.ID, layout=c(1, 5))


###################################################
### code chunk number 75: hab5b
###################################################
qqmath(~ resid(hd.lme.5) | stage$Hab.ID,
            prepanel = prepanel.qqmathline,
            panel = function(x, ...) {
               panel.qqmathline(x, distribution = qnorm)
               panel.qqmath(x, ...)
            })


###################################################
### code chunk number 76: habitat5a
###################################################
print(
plot(hd.lme.5, resid(.) ~ fitted(.) | Hab.ID, layout=c(1, 5))
      )


###################################################
### code chunk number 77: habitat5b
###################################################
print(
qqmath(~ resid(hd.lme.5) | stage$Hab.ID,
            prepanel = prepanel.qqmathline,
            panel = function(x, ...) {
               panel.qqmathline(x, distribution = qnorm)
               panel.qqmath(x, ...)
            })
)


###################################################
### code chunk number 78: hd.tree.5
###################################################
plot(augPred(hd.lme.5), 
     index.cond = list(panel.order),
     strip = strip.custom(par.strip.text = list(cex = 0.5)))


###################################################
### code chunk number 79: hd-tree-5
###################################################
print(
plot(augPred(hd.lme.5), 
     index.cond = list(panel.order),
     strip = strip.custom(par.strip.text = list(cex = 0.5)))
      )


###################################################
### code chunk number 80: addedv
###################################################
age.lme.1 <- lme(Age ~ dbhib.cm, 
                 random = ~1 | Forest.ID/Tree.ID,
                 data = stage)
res.Age <- residuals(age.lme.1, level = 0)
res.HD <- residuals(hd.lme.5, level = 0)
scatter.smooth(res.Age, res.HD, 
    xlab = "Variation unique to Age",
    ylab = "Variation in Height after all but Age")


###################################################
### code chunk number 81: added
###################################################
age.lme.1 <- lme(Age ~ dbhib.cm, 
                 random = ~1 | Forest.ID/Tree.ID,
                 data = stage)
res.Age <- residuals(age.lme.1, level = 0)
res.HD <- residuals(hd.lme.5, level = 0)
scatter.smooth(res.Age, res.HD, 
    xlab = "Variation unique to Age",
    ylab = "Variation in Height after all but Age")


###################################################
### code chunk number 82: addedhab
###################################################
xyplot(stage$height.m ~ fitted(hd.lme.5, level=0) | Hab.ID,
        xlab="Predicted height (m)",
        ylab="Observed height (m)",
        data=stage,
        panel = function(x, y, subscripts) {
                panel.xyplot(x, y)
                panel.abline(0, 1)
                panel.abline(lm(y ~ x), lty=3)
        }
)


###################################################
### code chunk number 83: hd-hab
###################################################
print(
xyplot(stage$height.m ~ fitted(hd.lme.5, level=0) | Hab.ID,
        xlab="Predicted height (m)",
        ylab="Observed height (m)",
        data=stage,
        panel = function(x, y, subscripts) {
                panel.xyplot(x, y)
                panel.abline(0, 1)
                panel.abline(lm(y ~ x), lty=3)
        }
)
      )


###################################################
### code chunk number 84: get.gutten
###################################################
Stangle("../../ch2/sweave/gutten.rnw")
source("gutten.R")


###################################################
### code chunk number 85: fitting_models.rnw:2422-2423
###################################################
library(lattice)


###################################################
### code chunk number 86: fig-gutten-data
###################################################
xyplot(dbh.cm ~ age.bh | tree.ID, type="l", data=gutten)


###################################################
### code chunk number 87: gutten-data
###################################################
print(
xyplot(dbh.cm ~ age.bh | tree.ID, type="l", data=gutten)
)


###################################################
### code chunk number 88: fitting_models.rnw:2442-2445
###################################################
library(nlme)
gutten.d <- groupedData(dbh.cm ~ age.bh | tree.ID, 
                        data = gutten)


###################################################
### code chunk number 89: fitting_models.rnw:2450-2453
###################################################
gutten.nlsList <- 
  nlsList(dbh.cm ~ SSasympOrig(age.bh, asymptote, scale), 
          data = gutten.d)


###################################################
### code chunk number 90: fig-gutten-residuals
###################################################
plot(gutten.nlsList, 
     residuals(., type="pearson") ~ fitted(.) | tree.ID)


###################################################
### code chunk number 91: gutten-residuals
###################################################
print(
plot(gutten.nlsList, 
     residuals(., type="pearson") ~ fitted(.) | tree.ID)
)


###################################################
### code chunk number 92: fig-gutten-intervals
###################################################
plot(intervals(gutten.nlsList), layout=c(2,1))


###################################################
### code chunk number 93: gutten-intervals
###################################################
print(
plot(intervals(gutten.nlsList), layout=c(2,1))
)


###################################################
### code chunk number 94: fitting_models.rnw:2509-2510
###################################################
methods(class=class(gutten.nlsList))


###################################################
### code chunk number 95: fitting_models.rnw:2515-2516
###################################################
methods(class=class(summary(gutten.nlsList)))


###################################################
### code chunk number 96: fitting_models.rnw:2520-2521
###################################################
str(coef(summary(gutten.nlsList)))


###################################################
### code chunk number 97: fitting_models.rnw:2525-2529
###################################################
asymptote <- 
  coef(summary(gutten.nlsList))[,"Estimate","asymptote"]
half.age <- log(2) /
  exp(coef(summary(gutten.nlsList))[,"Estimate","scale"])


###################################################
### code chunk number 98: fig-nlsList-plot
###################################################
opar <- par(las=1, mar=c(4,4,1,1))
scatter.smooth(half.age, asymptote, 
               xlab = "Age", ylab = "Asymptote")
par(opar)


###################################################
### code chunk number 99: nlsList-plot
###################################################
opar <- par(las=1, mar=c(4,4,1,1))
scatter.smooth(half.age, asymptote, 
               xlab = "Age", ylab = "Asymptote")
par(opar)


###################################################
### code chunk number 100: fitting_models.rnw:2600-2601
###################################################
gutten.nlme.0 <- nlme(gutten.nlsList)


###################################################
### code chunk number 101: fitting_models.rnw:2607-2613 (eval = FALSE)
###################################################
## gutten.nlme.0 <- 
##   nlme(dbh.cm ~ SSasympOrig(age.bh, asymptote, scale),
##        fixed = asymptote + scale ~ 1,
##        random = asymptote + scale ~ 1,
##        start = c(asymptote = 50, scale = -5),
##        data = gutten.d)


###################################################
### code chunk number 102: fig-nlme-ac-0
###################################################
plot(ACF(gutten.nlme.0, form = ~1|tree.ID), alpha=0.01)


###################################################
### code chunk number 103: fitting_models.rnw:2637-2639
###################################################
gutten.nlme.1 <- update(gutten.nlme.0, 
                         correlation = corARMA(p = 1, q = 2))


###################################################
### code chunk number 104: fig-nlme-ac-1
###################################################
plot(ACF(gutten.nlme.1, resType = "n", form = ~1|tree.ID), 
         alpha = 0.01)


###################################################
### code chunk number 105: nlme-ac-0
###################################################
print(
plot(ACF(gutten.nlme.0, form = ~1|tree.ID), alpha=0.01)
)


###################################################
### code chunk number 106: nlme-ac-1
###################################################
print(
plot(ACF(gutten.nlme.1, resType = "n", form = ~1|tree.ID), 
         alpha = 0.01)
)


###################################################
### code chunk number 107: fitting_models.rnw:2684-2690
###################################################
gutten.nlme.2 <- 
  nlme(dbh.cm ~ SSasymp(age.bh, asymptote, scale, R0),
       fixed = asymptote + scale + R0 ~ 1,
       random = asymptote + scale  + R0 ~ 1,
       start = c(asymptote = 40, scale = -0.07, R0 = -4),
       data = gutten.d)


###################################################
### code chunk number 108: fig-gutten-nlme
###################################################
plot(comparePred(gutten.nlme.0, gutten.nlme.2))


###################################################
### code chunk number 109: gutten-nlme
###################################################
print(
plot(comparePred(gutten.nlme.0, gutten.nlme.2))
)


###################################################
### code chunk number 110: fitting_models.rnw:2753-2754
###################################################
intervals(gutten.nlme.1)


###################################################
### code chunk number 111: fitting_models.rnw:2780-2786 (eval = FALSE)
###################################################
## (gutten.gamm <- gamm(dbh.cm ~ s(age.bh), 
##                     random=list(tree.ID = ~1),
##                     data=gutten.d))
## 
## plot(gutten.gamm$lme, residuals(., type="pearson") ~ fitted(.) | tree.ID)
## 


###################################################
### code chunk number 112: gutten.nlme.2
###################################################
try(
gutten.nlme.2 <- update(gutten.nlme.0,
                        fixed = list(asymptote ~ 1, scale ~ site),
                        start = c(asymptote=50, scale=rep(-5,5)))
)


###################################################
### code chunk number 113: fitting_models.rnw:2916-2918
###################################################
system("rm -fr package-Ch6")
package.skeleton(name = "package-Ch6")


