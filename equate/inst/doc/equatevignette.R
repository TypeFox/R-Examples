### R code from vignette source 'equatevignette.Rnw'

###################################################
### code chunk number 1: equatevignette.Rnw:375-377
###################################################
op <- options()
options(prompt = "R>", show.signif.stars = FALSE, warn = -1)


###################################################
### code chunk number 2: equatevignette.Rnw:379-382
###################################################
library("equate")
act.x <- as.freqtab(ACTmath[, 1:2])
act.y <- as.freqtab(ACTmath[, c(1, 3)])


###################################################
### code chunk number 3: equatevignette.Rnw:385-387
###################################################
head(act.x)
rbind(x = summary(act.x), y = summary(act.y))


###################################################
### code chunk number 4: equatevignette.Rnw:390-392
###################################################
neat.x <- freqtab(KBneat$x, scales = list(0:36, 0:12))
neat.y <- freqtab(KBneat$y, scales = list(0:36, 0:12))


###################################################
### code chunk number 5: equatevignette.Rnw:396-406
###################################################
attach(PISA)
r3items <- paste(items$itemid[items$clusterid == "r3a"])
r6items <- paste(items$itemid[items$clusterid == "r6"])
r5items <- paste(items$itemid[items$clusterid == "r5"])
r7items <- paste(items$itemid[items$clusterid == "r7"])
pisa <- freqtab(students[students$book == 6, ],
	items = list(c(r3items, r6items), c(r5items, r7items)),
	scales = list(0:31, 0:29))
round(data.frame(summary(pisa),
	row.names = c("r3r6", "r5r7")), 2)


###################################################
### code chunk number 6: equatevignette.Rnw:410-412
###################################################
plot(x = act.x, lwd = 2, xlab = "Score", ylab = "Count")
plot(neat.x)


###################################################
### code chunk number 7: plotunivar
###################################################
plot(x = act.x, lwd = 2, xlab = "Score", ylab = "Count")


###################################################
### code chunk number 8: plotbivar
###################################################
plot(neat.x)


###################################################
### code chunk number 9: equatevignette.Rnw:422-423
###################################################
plot(x = act.x, lwd = 2, xlab = "Score", ylab = "Count")


###################################################
### code chunk number 10: equatevignette.Rnw:432-433
###################################################
plot(neat.x)


###################################################
### code chunk number 11: equatevignette.Rnw:444-446 (eval = FALSE)
###################################################
## presmoothing(~ poly(total, 3, raw = T) + poly(anchor, 3, raw = T) +
## 	total:anchor, data = neat.x)


###################################################
### code chunk number 12: equatevignette.Rnw:451-454 (eval = FALSE)
###################################################
## neat.xsf <- with(as.data.frame(neat.x), cbind(total, total^2,
## 	total^3, anchor, anchor^2, anchor^3, total*anchor))
## presmoothing(neat.x, smooth = "loglinear", scorefun = neat.xsf)


###################################################
### code chunk number 13: equatevignette.Rnw:459-460
###################################################
neat.xs <- presmoothing(neat.x, smooth = "log", degrees = list(3, 1))


###################################################
### code chunk number 14: equatevignette.Rnw:465-470
###################################################
neat.xsmat <- presmoothing(neat.x, smooth = "loglinear",
	degrees = list(3, 1), stepup = TRUE)
plot(neat.xs)
plot(neat.x, neat.xsmat, ylty = 1:4)
round(rbind(x = summary(neat.x), xs = summary(neat.xs)), 2)


###################################################
### code chunk number 15: plotbivarsmooth1
###################################################
plot(neat.xs)


###################################################
### code chunk number 16: plotbivarsmooth2
###################################################
plot(neat.x, neat.xsmat, ylty = 1:5)


###################################################
### code chunk number 17: equatevignette.Rnw:480-481
###################################################
plot(neat.xs)


###################################################
### code chunk number 18: equatevignette.Rnw:490-491
###################################################
plot(neat.x, neat.xsmat, ylty = 1:5)


###################################################
### code chunk number 19: equatevignette.Rnw:501-503
###################################################
presmoothing(neat.x, smooth = "loglinear",
	degrees = list(c(3, 3), c(1, 1)), compare = TRUE)


###################################################
### code chunk number 20: equatevignette.Rnw:508-509
###################################################
equate(act.x, act.y, type = "mean")


###################################################
### code chunk number 21: equatevignette.Rnw:512-514
###################################################
neat.ef <- equate(neat.x, neat.y, type = "equip",
	method = "frequency estimation", smoothmethod = "log")


###################################################
### code chunk number 22: equatevignette.Rnw:519-520
###################################################
summary(neat.ef)


###################################################
### code chunk number 23: equatevignette.Rnw:524-526
###################################################
cbind(newx = c(3, 29, 8, 7, 13),
	yx = equate(c(3, 29, 8, 7, 13), y = neat.ef))


###################################################
### code chunk number 24: equatevignette.Rnw:529-530
###################################################
head(neat.ef$concordance)


###################################################
### code chunk number 25: equatevignette.Rnw:534-540
###################################################
neat.i <- equate(neat.x, neat.y, type = "ident")
neat.lt <- equate(neat.x, neat.y, type = "linear",
	method = "tucker")
neat.comp <- composite(list(neat.i, neat.lt), wc = .5,
	symmetric = TRUE)
plot(neat.comp, addident = FALSE)


###################################################
### code chunk number 26: plotcomposite
###################################################
plot(neat.comp, addident = FALSE)


###################################################
### code chunk number 27: equatevignette.Rnw:549-550
###################################################
plot(neat.comp, addident = FALSE)


###################################################
### code chunk number 28: equatevignette.Rnw:563-571
###################################################
pisa.i <- equate(pisa, type = "ident", lowp = c(3.5, 2))
pisa.m <- equate(pisa, type = "mean", lowp = c(3.5, 2))
pisa.l <- equate(pisa, type = "linear", lowp = c(3.5, 2))
pisa.c <- equate(pisa, type = "circ", lowp = c(3.5, 2))
pisa.e <- equate(pisa, type = "equip", smooth = "log",
	lowp = c(3.5, 2))
plot(pisa.i, pisa.m, pisa.l, pisa.c, pisa.e, addident = F,
	xpoints = pisa, morepars = list(ylim = c(0, 31)))


###################################################
### code chunk number 29: plotstudy2
###################################################
plot(pisa.i, pisa.m, pisa.l, pisa.c, pisa.e, addident = F,
	xpoints = pisa, morepars = list(ylim = c(0, 31)))


###################################################
### code chunk number 30: equatevignette.Rnw:579-580
###################################################
plot(pisa.i, pisa.m, pisa.l, pisa.c, pisa.e, addident = F,
	xpoints = pisa, morepars = list(ylim = c(0, 31)))


###################################################
### code chunk number 31: equatevignette.Rnw:591-595
###################################################
pisa.x <- freqtab(totals$b4[1:200, c("r3a", "r2", "s2")],
	scales = list(0:15, 0:17, 0:18))
pisa.y <- freqtab(totals$b4[201:400, c("r4a", "r2", "s2")],
	scales = list(0:16, 0:17, 0:18))


###################################################
### code chunk number 32: equatevignette.Rnw:600-606
###################################################
pisa.mnom <- equate(pisa.x, pisa.y, type = "mean",
	method = "nom")
pisa.mtuck <- equate(pisa.x, pisa.y, type = "linear",
	method = "tuck")
pisa.mfreq <- equate(pisa.x, pisa.y, type = "equip",
	method = "freq", smooth = "loglin")


###################################################
### code chunk number 33: equatevignette.Rnw:609-615
###################################################
pisa.snom <- equate(margin(pisa.x, 1:2), margin(pisa.y, 1:2),
	type = "mean", method = "nom")
pisa.stuck <- equate(margin(pisa.x, 1:2), margin(pisa.y, 1:2),
	type = "linear", method = "tuck")
pisa.sfreq <- equate(margin(pisa.x, 1:2), margin(pisa.y, 1:2),
	type = "equip", method = "freq", smooth = "loglin")


###################################################
### code chunk number 34: plotstudy3
###################################################
plot(pisa.snom, pisa.stuck, pisa.sfreq,
	pisa.mnom, pisa.mtuck, pisa.mfreq,
	col = rep(rainbow(3), 2), lty = rep(1:2, each = 3))


###################################################
### code chunk number 35: equatevignette.Rnw:625-626
###################################################
plot(pisa.snom, pisa.stuck, pisa.sfreq,
	pisa.mnom, pisa.mtuck, pisa.mfreq,
	col = rep(rainbow(3), 2), lty = rep(1:2, each = 3))


###################################################
### code chunk number 36: equatevignette.Rnw:643-651
###################################################
neat.xp <- presmoothing(neat.x, "loglinear", degrees = list(4, 2))
neat.xpmat <- presmoothing(neat.x, "loglinear", degrees = list(4, 2),
	stepup = TRUE)
neat.yp <- presmoothing(neat.y, "loglinear", degrees = list(4, 2))
neat.ypmat <- presmoothing(neat.y, "loglinear", degrees = list(4, 2),
	stepup = TRUE)
plot(neat.x, neat.xpmat)
plot(neat.y, neat.ypmat)


###################################################
### code chunk number 37: plotstudy1x
###################################################
plot(neat.x, neat.xpmat)


###################################################
### code chunk number 38: plotstudy1y
###################################################
plot(neat.y, neat.ypmat)


###################################################
### code chunk number 39: equatevignette.Rnw:662-663
###################################################
plot(neat.x, neat.xpmat)


###################################################
### code chunk number 40: equatevignette.Rnw:672-673
###################################################
plot(neat.y, neat.ypmat)


###################################################
### code chunk number 41: equatevignette.Rnw:681-686
###################################################
set.seed(131031)
reps <- 100
xn <- 100
yn <- 100
crit <- equate(neat.xp, neat.yp, "e", "c")$conc$yx


###################################################
### code chunk number 42: equatevignette.Rnw:689-700
###################################################
neat.args <- list(i = list(type = "i"),
	mt = list(type = "mean", method = "t"),
	mc = list(type = "mean", method = "c"),
	lt = list(type = "lin", method = "t"),
	lc = list(type = "lin", method = "c"),
	ef = list(type = "equip", method = "f", smooth = "log"),
	ec = list(type = "equip", method = "c", smooth = "log"),
	ct = list(type = "circ", method = "t"),
	cc = list(type = "circ", method = "c", chainmidp = "lin"))
bootout <- bootstrap(x = neat.xp, y = neat.yp, xn = xn, yn = yn,
	reps = reps, crit = crit, args = neat.args)


###################################################
### code chunk number 43: equatevignette.Rnw:703-710
###################################################
plot(bootout, addident = F, col = c(1, rainbow(8)))
plot(bootout, out = "se", addident = F,
	col = c(1, rainbow(8)), legendplace = "top")
plot(bootout, out = "bias", addident = F, legendplace = "top",
	col = c(1, rainbow(8)), morepars = list(ylim = c(-.9, 3)))
plot(bootout, out = "rmse", addident = F, legendplace = "top",
	col = c(1, rainbow(8)), morepars = list(ylim = c(0, 3)))


###################################################
### code chunk number 44: plotstudy1means
###################################################
plot(bootout, addident = F, col = c(1, rainbow(8)))


###################################################
### code chunk number 45: plotstudy1se
###################################################
plot(bootout, out = "se", addident = F,
	col = c(1, rainbow(8)), legendplace = "top")


###################################################
### code chunk number 46: plotstudy1bias
###################################################
plot(bootout, out = "bias", addident = F, legendplace = "top",
	col = c(1, rainbow(8)), morepars = list(ylim = c(-.9, 3)))


###################################################
### code chunk number 47: plotstudy1rmse
###################################################
plot(bootout, out = "rmse", addident = F, legendplace = "top",
	col = c(1, rainbow(8)), morepars = list(ylim = c(0, 3)))


###################################################
### code chunk number 48: equatevignette.Rnw:729-730
###################################################
plot(bootout, addident = F, col = c(1, rainbow(8)))


###################################################
### code chunk number 49: equatevignette.Rnw:739-740
###################################################
plot(bootout, out = "se", addident = F,
	col = c(1, rainbow(8)), legendplace = "top")


###################################################
### code chunk number 50: equatevignette.Rnw:749-750
###################################################
plot(bootout, out = "bias", addident = F, legendplace = "top",
	col = c(1, rainbow(8)), morepars = list(ylim = c(-.9, 3)))


###################################################
### code chunk number 51: equatevignette.Rnw:759-760
###################################################
plot(bootout, out = "rmse", addident = F, legendplace = "top",
	col = c(1, rainbow(8)), morepars = list(ylim = c(0, 3)))


###################################################
### code chunk number 52: equatevignette.Rnw:768-769
###################################################
round(summary(bootout), 2)


###################################################
### code chunk number 53: equatevignette.Rnw:771-773
###################################################
detach(PISA)
options(op)


