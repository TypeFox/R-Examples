### R code from vignette source 'etmCIF_tutorial.Rnw'

###################################################
### code chunk number 1: etmCIF_tutorial.Rnw:34-36
###################################################
require(etm)
data(abortion)


###################################################
### code chunk number 2: etmCIF_tutorial.Rnw:50-51
###################################################
head(abortion)


###################################################
### code chunk number 3: etmCIF_tutorial.Rnw:95-98
###################################################
cif.abortion <- etmCIF(Surv(entry, exit, cause != 0) ~ group, 
                   abortion, etype = cause, failcode = 3)
cif.abortion


###################################################
### code chunk number 4: etmCIF_tutorial.Rnw:107-108
###################################################
s.cif.ab <- summary(cif.abortion)


###################################################
### code chunk number 5: etmCIF_tutorial.Rnw:115-116
###################################################
s.cif.ab


###################################################
### code chunk number 6: etmCIF_tutorial.Rnw:128-129
###################################################
plot(cif.abortion)


###################################################
### code chunk number 7: etmCIF_tutorial.Rnw:144-147
###################################################
plot(cif.abortion, curvlab = c("Control", "Exposed"), ylim = c(0, 0.6), 
     ci.type = "bars", pos.ci = 27, col = c(1, 2), ci.lwd = 6, 
     lwd = 2, lty = 1, cex = 1.3)


###################################################
### code chunk number 8: etmCIF_tutorial.Rnw:166-169
###################################################
plot(cif.abortion, curvlab = c("Control", "Exposed"), ylim = c(0, 0.6), 
     ci.type = "bars", pos.ci = c(27, 28), col = c(1, 1), ci.lwd = 6, 
     lwd = 2, lty = c(2, 1), cex = 1.3)


###################################################
### code chunk number 9: etmCIF_tutorial.Rnw:182-184
###################################################
plot(cif.abortion, curvlab = c("Control", "Exposed"), ylim = c(0, 0.5), 
     ci.type = "pointwise", col = c(1, 2), lwd = 2, lty = 1, cex = 1.3)


###################################################
### code chunk number 10: etmCIF_tutorial.Rnw:199-205
###################################################
plot(cif.abortion, which.cif = c(1, 2), ylim = c(0, 0.8), lwd = 2,
     col = c(1, 1, 2, 2), lty = c(1, 2, 1, 2), legend = FALSE)
legend(0, 0.8, c("Control", "Exposed"), col = c(1, 2), lty = 1, 
       bty = "n", lwd = 2)
legend(0, 0.7, c("ETOP", "Life Birth"), col = 1, lty = c(1, 2), 
       bty = "n", lwd = 2)


###################################################
### code chunk number 11: etmCIF_tutorial.Rnw:225-231
###################################################
abortion$status <- with(abortion, ifelse(cause == 2, "life birth",
                        ifelse(cause == 1, "ETOP", "spontaneous abortion")))
abortion$status <- factor(abortion$status)

abortion$treat <- with(abortion, ifelse(group == 0, "control", "exposed"))
abortion$treat <- factor(abortion$treat)


###################################################
### code chunk number 12: etmCIF_tutorial.Rnw:236-239
###################################################
new.cif <- etmCIF(Surv(entry, exit, status != 0) ~ treat, abortion, 
                  etype = status, failcode = "spontaneous abortion")
new.cif


###################################################
### code chunk number 13: etmCIF_tutorial.Rnw:260-261
###################################################
trprob(new.cif[[1]], "0 spontaneous abortion", c(1, 10, 27))


###################################################
### code chunk number 14: etmCIF_tutorial.Rnw:275-276 (eval = FALSE)
###################################################
## lines(cif.abortion[[2]], tr.choice = "0 1", col = 2, lwd = 2)


###################################################
### code chunk number 15: etmCIF_tutorial.Rnw:281-285
###################################################
plot(cif.abortion, curvlab = c("Control", "Exposed"), ylim = c(0, 0.6), 
     ci.type = "bars", pos.ci = c(27, 28), col = c(1, 1), ci.lwd = 6, 
     lwd = 2, lty = c(2, 1), cex = 1.3)
lines(cif.abortion[[2]], tr.choice = "0 1", col = 2, lwd = 2)


