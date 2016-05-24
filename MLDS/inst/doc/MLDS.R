### R code from vignette source 'MLDS.Stex'

###################################################
### code chunk number 1: MLDS.Stex:72-73
###################################################
library(MLDS)


###################################################
### code chunk number 2: MLDS.Stex:252-253
###################################################
head(kk3)


###################################################
### code chunk number 3: MLDS.Stex:699-700
###################################################
c(resp = 1, S1 = 7, S2 = 9, S3 = 2, S4 = 4)


###################################################
### code chunk number 4: MLDS.Stex:795-799
###################################################
data(AutumnLab)
x.mlds <- mlds(AutumnLab)
x.6pt <- Get6pts(x.mlds, nrep = 1)
t(sapply(x.6pt, function(x) x[4, ]))


###################################################
### code chunk number 5: MLDS.Stex:806-807
###################################################
unlist(attr(x.6pt, 'indices')[4, ])


###################################################
### code chunk number 6: MLDS.Stex:853-857
###################################################
data(kk1)
data(kk2)
data(kk3)
kk <- SwapOrder(rbind(kk1, kk2, kk3))


###################################################
### code chunk number 7: MLDS.Stex:864-868
###################################################
kk.mlds <- mlds(kk)
summary(kk.mlds)
kkopt.mlds <- mlds(kk, method = "optim", opt.init = c(seq(0, 1, len = 11), 0.2))
summary(kkopt.mlds)


###################################################
### code chunk number 8: MLDS.Stex:998-1000
###################################################
kk.mlds.logit <- mlds(kk, lnk = "logit")
kk.mlds.cauchit <- mlds(kk, lnk = "cauchit")


###################################################
### code chunk number 9: MLDS.Stex:1104-1105
###################################################
kk[residuals(kk.mlds$obj) < -2.5, ]


###################################################
### code chunk number 10: MLDS.Stex:1137-1139 (eval = FALSE)
###################################################
## library(brglm)
## mlds(kk2, glm.meth = brglm.fit)


###################################################
### code chunk number 11: MLDS.Stex:1187-1188
###################################################
kk.frm <- mlds(~ (sx/0.98)^p[1], p = c(2, 0.2), data = kk)


###################################################
### code chunk number 12: MLDS.Stex:1212-1214
###################################################
c(p = kk.frm$par, sigma = kk.frm$sigma)
sqrt(diag(solve(kk.frm$hess)))


###################################################
### code chunk number 13: MLDS.Stex:1226-1230
###################################################
ddf <- diff(c(attr(logLik(kk.frm), "df"), 
	attr(logLik(kk.mlds), "df")))
pchisq(-2 * c(logLik(kk.frm) - logLik(kk.mlds)), 
	ddf, lower.tail = FALSE)


###################################################
### code chunk number 14: MLDS.Stex:1246-1251
###################################################
plot(kk.mlds, standard.scale = TRUE,
	xlab = expression(r^2),
	ylab = "Difference Scale Value")
xx <- seq(0, 0.98, len = 100)
lines(xx, kk.frm$func(kk.frm$par, xx))


