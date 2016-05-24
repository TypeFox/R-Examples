### R code from vignette source 'ExtExamp/ExtExamp.Stex'

###################################################
### code chunk number 1: ExtExamp.Stex:9-28
###################################################
options(width=60, digits = 3, str = strOptions(strict.width = "cut"))
pkg <- search()[2] 
while (search()[2] !=  if(.Platform$GUI == "AQUA") "tools:RGUI" else "package:stats") {
#detach(pos = match(pkg, search()))
spkg <- strsplit( pkg, ":"  )[[1]][2] 
if (packageHasNamespace(spkg, .libPaths()[1]) )
unloadNamespace(spkg ) else detach(pos = match(pkg, search()))
 pkg <- search()[2]
 }
rm(list = ls())
gc()
library(MPDiR)
#library(lattice) 
#ltheme <- canonical.theme("pdf", color = FALSE) ## in-built B&W theme 
#ltheme$strip.background$bg <- "grey85" ## change strip bg 
#lattice.options(default.theme = ltheme) ## set as default 
formals(deparse)$width.cutoff <- 50L
assignInNamespace("deparse", deparse, "base")
rm(deparse)


###################################################
### code chunk number 2: ExtExamp.Stex:90-91 (eval = FALSE)
###################################################
## install.packages("MPDiR")


###################################################
### code chunk number 3: ExtExamp.Stex:97-98
###################################################
library(MPDiR)


###################################################
### code chunk number 4: ExtExamp.Stex:112-113
###################################################
 search()


###################################################
### code chunk number 5: ExtExamp.Stex:123-124
###################################################
data(HSP)


###################################################
### code chunk number 6: ExtExamp.Stex:128-129
###################################################
ls()


###################################################
### code chunk number 7: ExtExamp.Stex:141-142
###################################################
str(HSP)


###################################################
### code chunk number 8: ExtExamp.Stex:176-177
###################################################
names(HSP)


###################################################
### code chunk number 9: ExtExamp.Stex:183-184
###################################################
names(HSP) <- c("Quanta", "PerCent", "N", "Obs", "Run")


###################################################
### code chunk number 10: ExtExamp.Stex:197-198
###################################################
head(HSP, n = 10)


###################################################
### code chunk number 11: ExtExamp.Stex:211-212
###################################################
summary(HSP)


###################################################
### code chunk number 12: ExtExamp.Stex:220-221
###################################################
with(HSP, table(Run, Obs))


###################################################
### code chunk number 13: ExtExamp.Stex:241-242
###################################################
   plot(HSP)


###################################################
### code chunk number 14: ExtExamp.Stex:276-277 (eval = FALSE)
###################################################
## methods(plot)


###################################################
### code chunk number 15: ExtExamp.Stex:299-306 (eval = FALSE)
###################################################
## xyplot(PerCent/100 ~ Quanta | Obs + Run, data = HSP, 
##   xlab = "Quanta/flash", ylab = "Proportion \"seen\"",
##   scales = list(x = list(log = TRUE, limits = c(10, 500), 
##     at = c(10, 20, 50, 100, 200),
##     labels =  c(10, 20, 50, 100, 200))),
##   skip = c(rep(F, 5), T), layout = c(3, 2, 1),
##   as.table = TRUE)


###################################################
### code chunk number 16: ExtExamp.Stex:309-315
###################################################
library(lattice) 
ltheme <- canonical.theme(color = FALSE) ## in-built B&W theme 
ltheme$strip.background$bg <- "grey85" ## change strip bg 
ltheme$strip.background$col <- "grey85" ## change strip bg 
lattice.options(default.theme = ltheme) ## set as default 
rm(ltheme)


###################################################
### code chunk number 17: ExtExamp.Stex:319-332
###################################################
library(lattice) 
ltheme <- canonical.theme("pdf", color = FALSE) ## in-built B&W theme 
ltheme$strip.background$bg <- "grey85" ## change strip bg 
lattice.options(default.theme = ltheme) ## set as default 
rm(ltheme)
print(
xyplot(PerCent/100  ~ Quanta | Obs + Run, data = HSP, 
	xlab = "Quanta/flash", ylab = "Proportion \"seen\"",
	scales = list(x = list(log = TRUE, limits = c(10, 500), at = c(10, 20, 50, 100, 200),
		labels =  c(10, 20, 50, 100, 200))),
	skip = c(rep(F, 5), T), layout = c(3, 2, 1),
	as.table = TRUE)
)


###################################################
### code chunk number 18: ExtExamp.Stex:451-453
###################################################
HSP$NumYes <- round(HSP$N * HSP$PerCent/100)
HSP$NumNo <- HSP$N - HSP$NumYes


###################################################
### code chunk number 19: ExtExamp.Stex:464-466
###################################################
HSP <- within(HSP, {NumYes <- round(N * PerCent/100)
    NumNo <- N - NumYes })


###################################################
### code chunk number 20: ExtExamp.Stex:482-483
###################################################
( SHR1 <- subset(HSP, Obs == "SH" & Run == "R1") )


###################################################
### code chunk number 21: ExtExamp.Stex:495-496 (eval = FALSE)
###################################################
## with(HSP, HSP[Obs == "SH" & Run == "R1", ])


###################################################
### code chunk number 22: ExtExamp.Stex:508-510
###################################################
SHR1.glm <- glm(formula = cbind(NumYes, NumNo) ~ log(Quanta), 
	family = binomial(probit), data = SHR1)


###################################################
### code chunk number 23: ExtExamp.Stex:545-546
###################################################
 SHR1.glm


###################################################
### code chunk number 24: ExtExamp.Stex:587-588
###################################################
anova(SHR1.glm, test = "Chisq")


###################################################
### code chunk number 25: ExtExamp.Stex:600-601
###################################################
summary(SHR1.glm)


###################################################
### code chunk number 26: ExtExamp.Stex:647-649
###################################################
thresh <- exp(qnorm(p = 0.6, mean = -coef(SHR1.glm)[1]/coef(SHR1.glm)[2], 
	sd = 1/coef(SHR1.glm)[2]))


###################################################
### code chunk number 27: ExtExamp.Stex:651-665 (eval = FALSE)
###################################################
## plot(PerCent/100 ~ Quanta, SHR1, 
##   xlab = "Quanta/Flash", ylab = "Proportion \"seen\"",  
##   main = "Obs: SH, Run:  1", log = "x",
##   xlim = c(20, 440))
## xseq <- seq(20, 450, len = 100)
## SHR1.pred <- predict(SHR1.glm, newdata = 
##    data.frame(Quanta = xseq), type = "response",
##    se.fit = TRUE)
## polygon(c(xseq, rev(xseq)), c(SHR1.pred$fit + SHR1.pred$se.fit, 
## 	rev(SHR1.pred$fit - SHR1.pred$se.fit)), 
## 	border = "white", col = "grey")
## lines(xseq, SHR1.pred$fit, lwd = 2)
## points(PerCent/100 ~ Quanta, SHR1, cex = 1.5, pch = 21, 
##   bg = "white")


###################################################
### code chunk number 28: ExtExamp.Stex:670-685
###################################################
 opar <- par(mfrow = c(1, 1))
plot(PerCent/100 ~ Quanta, SHR1, xlab = "Quant/Flash", ylab = "Proportion \"seen\"",
	main = "Obs: SH, Run:  1", log = "x", #cex = 1.5, 
	xlim = c(20, 440))
xseq <- seq(20, 450, len = 100)
SHR1.pred <- predict(SHR1.glm, newdata = data.frame(Quanta = xseq), type = "response",
   se.fit = TRUE)
polygon(c(xseq, rev(xseq)), c(SHR1.pred$fit + SHR1.pred$se.fit, 
	rev(SHR1.pred$fit - SHR1.pred$se.fit)), 
	border = "white", col = "grey")
lines(xseq, SHR1.pred$fit, lwd = 2)
points(PerCent/100 ~ Quanta, SHR1, cex = 1.5, pch = 21, 
  bg = "white")
segments(c(20, thresh), c(0.6, 0.6), c(thresh, thresh), c(0.6, 0), lty = 3, lwd = 3)
par(opar)


###################################################
### code chunk number 29: ExtExamp.Stex:699-700
###################################################
confint(SHR1.glm)


###################################################
### code chunk number 30: ExtExamp.Stex:723-726
###################################################
(thresh <- exp(qnorm(p = 0.6, 
	mean = -coef(SHR1.glm)[1]/coef(SHR1.glm)[2], 
	sd = 1/coef(SHR1.glm)[2])))


###################################################
### code chunk number 31: ExtExamp.Stex:743-749
###################################################
	thresh.est <- function(p, obj) {
	cc <- coef(obj)
	m <- -cc[1]/cc[2]
	std <- 1/cc[2]
	qnorm(p, m, std)
	}


###################################################
### code chunk number 32: ExtExamp.Stex:757-758
###################################################
thresh.est(0.6, SHR1.glm)


###################################################
### code chunk number 33: ExtExamp.Stex:761-762
###################################################
exp(thresh.est(0.6, SHR1.glm))


###################################################
### code chunk number 34: ExtExamp.Stex:770-771
###################################################
diff(thresh.est(c(0.5, 0.75), SHR1.glm))


###################################################
### code chunk number 35: ExtExamp.Stex:784-785
###################################################
HSP$id <- with(HSP, interaction(Obs, Run, drop = TRUE))


###################################################
### code chunk number 36: ExtExamp.Stex:795-796
###################################################
HSP0.glm <- glm(cbind(NumYes, NumNo) ~ log(Quanta), binomial(probit), HSP)


###################################################
### code chunk number 37: ExtExamp.Stex:802-803
###################################################
HSP1.glm <- update(HSP0.glm,  . ~ id/log(Quanta) - 1)


###################################################
### code chunk number 38: ExtExamp.Stex:820-821
###################################################
anova(HSP0.glm, HSP1.glm, test = "Chisq")


###################################################
### code chunk number 39: ExtExamp.Stex:835-857 (eval = FALSE)
###################################################
## nd <- data.frame(
##   Quanta = rep(seq(20, 500, len = 100), 5),
##   id = rep(levels(HSP$id), each = 100)
##   )
## levels(nd$id) <- levels(HSP$id)
## nd$pred0 <- predict(HSP0.glm, newdata = nd, type = "response")
## nd$pred1 <- predict(HSP1.glm, newdata = nd, type = "response")
## mm <- matrix(unlist(strsplit(as.character(nd$id), 
##   "\\.")), ncol = 2, byrow = TRUE)
## nd$Obs <- factor(mm[, 1], levels = c("SH", "SS", "MHP"))
## nd$Run <- factor(mm[, 2])
## 
## library(ggplot2)
## 
## qplot(Quanta, PerCent/100, data = HSP, facets = Run ~ Obs) + 
##   geom_point(size = 4) +
##   geom_line(data = nd, aes(x = Quanta, y = pred1), size = 1) +
##   geom_line(data = nd, aes(x = Quanta, y = pred0), size = 1,
## 	linetype = "dashed") +
##   scale_x_log10(limits = c(10, 500), 
##     breaks = c(10, 20, 50, 100, 200),
##     labels =  c(10, 20, 50, 100, 200))		


###################################################
### code chunk number 40: ExtExamp.Stex:861-884
###################################################
library(ggplot2)
nd <- data.frame(
	Quanta = rep(seq(20, 500, len = 100), 5),
	id = rep(levels(HSP$id), each = 100)
	)
levels(nd$id) <- levels(HSP$id)
nd$pred0 <- predict(HSP0.glm, newdata = nd, type = "response")
nd$pred1 <- predict(HSP1.glm, newdata = nd, type = "response")
mm <- matrix(unlist(strsplit(as.character(nd$id), 
			"\\.")), ncol = 2, byrow = TRUE)
nd$Obs <- factor(mm[, 1], levels = c("SH", "SS", "MHP"))
nd$Run <- factor(mm[, 2])

library(ggplot2)
print(
qplot(Quanta, PerCent/100, data = HSP, facets = Run ~ Obs) + 
geom_point(size = 4) +
geom_line(data = nd, aes(x = Quanta, y = pred1), size = 1) +
geom_line(data = nd, aes(x = Quanta, y = pred0), size = 1,
	linetype = "dashed") +
scale_x_log10(limits = c(10, 500), breaks = c(10, 20, 50, 100, 200),
	labels =  c(10, 20, 50, 100, 200))		
)


###################################################
### code chunk number 41: ExtExamp.Stex:913-914 (eval = FALSE)
###################################################
##  rm(list = ls())


###################################################
### code chunk number 42: ExtExamp.Stex:966-969 (eval = FALSE)
###################################################
## head(HSP, n = n)
## head(HSP, n = 15)
## head(HSP, n <- 15)


