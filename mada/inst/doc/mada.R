### R code from vignette source 'mada.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: mada.Rnw:94-95
###################################################
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)


###################################################
### code chunk number 2: mada.Rnw:107-108 (eval = FALSE)
###################################################
## install.packages("mada")


###################################################
### code chunk number 3: mada.Rnw:111-112
###################################################
library("mada")


###################################################
### code chunk number 4: mada.Rnw:144-146
###################################################
y <- 142 * .944 
y


###################################################
### code chunk number 5: mada.Rnw:149-150
###################################################
round(y)


###################################################
### code chunk number 6: mada.Rnw:155-160
###################################################
AuditC6 <- data.frame(TP = c(47, 126, 19, 36, 130, 84),
                      FN = c(9, 51, 10, 3, 19, 2),
                      FP = c(101, 272, 12, 78, 211, 68),
                      TN = c(738, 1543, 192, 276, 959, 89))
AuditC6


###################################################
### code chunk number 7: mada.Rnw:163-165
###################################################
AuditC6$names <- c("Study 1", "Study 2", "Study 4",
                   "Study 4", "Study 5", "Study 6")


###################################################
### code chunk number 8: mada.Rnw:169-171
###################################################
data("AuditC")
tail(AuditC)


###################################################
### code chunk number 9: mada.Rnw:192-193 (eval = FALSE)
###################################################
## madad(AuditC)


###################################################
### code chunk number 10: mada.Rnw:223-224 (eval = FALSE)
###################################################
## madad(AuditC, level = 0.80)


###################################################
### code chunk number 11: mada.Rnw:227-229
###################################################
AuditC.d <- madad(AuditC)
AuditC.d$fpr


###################################################
### code chunk number 12: mada.Rnw:245-247 (eval = FALSE)
###################################################
## forest(madad(AuditC), type = "sens")
## forest(madad(AuditC), type = "spec")


###################################################
### code chunk number 13: mada.Rnw:250-258
###################################################
pdf(file = "pairedforest.pdf", width = 12, height = 6)
par(omi = c(0,0,0,0), mai = c(0.9,0.3,0.3,0.3))
plot.new()
par(fig = c(0, 0.5, 0, 1), pty = "s", new = TRUE)
forest(AuditC.d, type = "sens", xlab = "Sensitivity")
par(fig = c(0.5, 1, 0, 1), pty = "s", new = TRUE)
forest(AuditC.d, type = "spec", xlab = "Specificity")
dev.off()


###################################################
### code chunk number 14: mada.Rnw:274-278
###################################################
rs <- rowSums(AuditC)
weights <- 4 * rs / max(rs)
crosshair(AuditC, xlim = c(0,0.6), ylim = c(0.4,1), 
          col = 1:14, lwd = weights)


###################################################
### code chunk number 15: mada.Rnw:282-284
###################################################
ROCellipse(AuditC, pch = "")
points(fpr(AuditC), sens(AuditC))


###################################################
### code chunk number 16: mada.Rnw:287-296
###################################################
pdf(file = "diagplots.pdf", width = 12, height = 6)
par(omi = c(0,0,0,0), mai = c(0.9,0.9,0.3,0.3))
plot.new()
par(fig = c(0, 0.5, 0, 1), pty = "s", new = TRUE)
crosshair(AuditC, xlim = c(0,0.6), ylim = c(0.4,1), col = 1:14, lwd = weights)
par(fig = c(0.5, 1, 0, 1), pty = "s", new = TRUE)
ROCellipse(AuditC, pch = "")
points(fpr(AuditC), sens(AuditC))
dev.off()


###################################################
### code chunk number 17: mada.Rnw:343-345
###################################################
(fit.DOR.DSL <- madauni(AuditC))
(fit.DOR.MH <- madauni(AuditC, method = "MH"))


###################################################
### code chunk number 18: mada.Rnw:348-349
###################################################
summary(fit.DOR.DSL)


###################################################
### code chunk number 19: mada.Rnw:352-353
###################################################
forest(fit.DOR.DSL)


###################################################
### code chunk number 20: mada.Rnw:356-360
###################################################
pdf(file = "DORforest.pdf", width = 6, height = 6)
par(omi = c(0,0,0,0), mai = c(0.9,0.9,0.3,0.3))
forest(fit.DOR.DSL)
dev.off()


###################################################
### code chunk number 21: mada.Rnw:377-379
###################################################
(fit.phm.homo <- phm(AuditC, hetero = FALSE))
(fit.phm.het <- phm(AuditC))


###################################################
### code chunk number 22: mada.Rnw:382-383
###################################################
summary(fit.phm.homo)


###################################################
### code chunk number 23: mada.Rnw:386-387
###################################################
summary(fit.phm.het)


###################################################
### code chunk number 24: mada.Rnw:390-392
###################################################
plot(fit.phm.het, xlim = c(0,0.6), ylim = c(0.4,1))
ROCellipse(AuditC, add = TRUE)


###################################################
### code chunk number 25: mada.Rnw:395-400
###################################################
pdf(file = "phmplot.pdf", width = 6, height = 6)
par(omi = c(0,0,0,0), mai = c(0.9,0.9,0.3,0.3))
plot(fit.phm.het, xlim = c(0,0.6), ylim = c(0.4,1))
ROCellipse(AuditC, add = TRUE)
dev.off()


###################################################
### code chunk number 26: mada.Rnw:446-447
###################################################
(fit.reitsma <- reitsma(AuditC))


###################################################
### code chunk number 27: mada.Rnw:450-451
###################################################
summary(fit.reitsma)


###################################################
### code chunk number 28: mada.Rnw:454-459
###################################################
plot(fit.reitsma, sroclwd = 2,
     main = "SROC curve (bivariate model) for AUDIT-C data")
points(fpr(AuditC), sens(AuditC), pch = 2)
legend("bottomright", c("data", "summary estimate"), pch = c(2,1))
legend("bottomleft", c("SROC", "conf. region"), lwd = c(2,1))


###################################################
### code chunk number 29: mada.Rnw:462-470
###################################################
pdf(file = "SROCAuditC.pdf", width = 6, height = 6)
par(omi = c(0,0,0,0), mai = c(0.9,0.9,0.3,0.3))
plot(fit.reitsma, sroclwd = 2,
     main = "SROC curve (bivariate model) for AUDIT-C data")
points(fpr(AuditC), sens(AuditC), pch = 2)
legend("bottomright", c("data", "summary estimate"), pch = c(2,1))
legend("bottomleft", c("SROC", "conf. region"), lwd = c(2,1))
dev.off()


###################################################
### code chunk number 30: mada.Rnw:484-488
###################################################
data("IAQ")
data("SAQ")
fit.IAQ <- reitsma(IAQ)
fit.SAQ <- reitsma(SAQ)


###################################################
### code chunk number 31: mada.Rnw:491-498
###################################################
plot(fit.IAQ, xlim = c(0,.5), ylim = c(.5,1),
     main = "Comparison of IAQ and SAQ")
lines(sroc(fit.SAQ), lty = 2)
ROCellipse(fit.SAQ, lty = 2, pch = 2, add = TRUE)
points(fpr(IAQ), sens(IAQ), cex = .5)
points(fpr(SAQ), sens(SAQ), pch = 2, cex = 0.5)
legend("bottomright", c("IAQ", "SAQ"), pch = 1:2, lty = 1:2)


###################################################
### code chunk number 32: mada.Rnw:500-510
###################################################
pdf(file = "SAQIAQ.pdf", width = 6, height = 6)
par(omi = c(0,0,0,0), mai = c(0.9,0.9,0.3,0.3))
plot(fit.IAQ, xlim = c(0,.5), ylim = c(.5,1),
     main = "Comparison of IAQ and SAQ")
lines(sroc(fit.SAQ), lty = 2)
ROCellipse(fit.SAQ, lty = 2, pch = 2, add = TRUE)
points(fpr(IAQ), sens(IAQ), cex = .5)
points(fpr(SAQ), sens(SAQ), pch = 2, cex = 0.5)
legend("bottomright", c("IAQ", "SAQ"), pch = 1:2, lty = 1:2)
dev.off()


###################################################
### code chunk number 33: mada.Rnw:523-524
###################################################
data("smoking")


###################################################
### code chunk number 34: mada.Rnw:527-528
###################################################
summary(smoking$type)


###################################################
### code chunk number 35: mada.Rnw:531-533
###################################################
fit.smoking.type <- reitsma(smoking, 
                            formula = cbind(tsens, tfpr) ~ type)


###################################################
### code chunk number 36: mada.Rnw:537-538
###################################################
summary(fit.smoking.type)


###################################################
### code chunk number 37: mada.Rnw:546-553
###################################################
fit.smoking.ml.type <- reitsma(smoking, 
                          formula = cbind(tsens, tfpr) ~ type, 
                          method = "ml")
fit.smoking.ml.intercept <- reitsma(smoking, 
                                    formula = cbind(tsens, tfpr) ~ 1,
                                    method = "ml")
anova(fit.smoking.ml.type, fit.smoking.ml.intercept)


###################################################
### code chunk number 38: mada.Rnw:563-569
###################################################
fit.smoking1 <- reitsma(smoking, method = "ml")
fit.smoking2 <- reitsma(smoking, 
                        alphasens = 0, alphafpr = 2, 
                        method = "ml")
AIC(fit.smoking1)
AIC(fit.smoking2)


