### R code from vignette source 'deTestSet.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
library("deTestSet")
options(prompt = "> ")
options(width=70)


###################################################
### code chunk number 2: deTestSet.Rnw:119-120
###################################################
out <- andrews()


###################################################
### code chunk number 3: andrews
###################################################
plot(out, which = 1:9, lwd = 2, ask = FALSE)
mtext(outer = TRUE, side = 3, line = -1.5, "andrews", cex = 1.5)


###################################################
### code chunk number 4: figandrews
###################################################
plot(out, which = 1:9, lwd = 2, ask = FALSE)
mtext(outer = TRUE, side = 3, line = -1.5, "andrews", cex = 1.5)


###################################################
### code chunk number 5: deTestSet.Rnw:142-143
###################################################
out <- beam()


###################################################
### code chunk number 6: beam
###################################################
plot(out, which = 1:9, lwd = 2, ask = FALSE)
mtext(outer = TRUE, side = 3, line = -1.5, "beam", cex = 1.5)


###################################################
### code chunk number 7: figbeam
###################################################
plot(out, which = 1:9, lwd = 2, ask = FALSE)
mtext(outer = TRUE, side = 3, line = -1.5, "beam", cex = 1.5)


###################################################
### code chunk number 8: deTestSet.Rnw:166-167
###################################################
out <- caraxis()


###################################################
### code chunk number 9: caraxis
###################################################
plot(out, which = 1:9, lwd = 2, ask = FALSE)
mtext(outer = TRUE, side = 3, line = -1.5, "caraxis", cex = 1.5)


###################################################
### code chunk number 10: figcaraxis
###################################################
plot(out, which = 1:9, lwd = 2, ask = FALSE)
mtext(outer = TRUE, side = 3, line = -1.5, "caraxis", cex = 1.5)


###################################################
### code chunk number 11: deTestSet.Rnw:190-191
###################################################
out <- crank()


###################################################
### code chunk number 12: crank
###################################################
plot(out, which = 1:9, lwd = 2, ask = FALSE)
mtext(outer = TRUE, side = 3, line = -1.5, "crank", cex = 1.5)


###################################################
### code chunk number 13: figcrank
###################################################
plot(out, which = 1:9, lwd = 2, ask = FALSE)
mtext(outer = TRUE, side = 3, line = -1.5, "crank", cex = 1.5)


###################################################
### code chunk number 14: deTestSet.Rnw:214-215
###################################################
out <- E5()


###################################################
### code chunk number 15: E5
###################################################
plot(out, lwd = 2, ask = FALSE)
mtext(outer = TRUE, side = 3, line = -1.5, "E5", cex = 1.5)


###################################################
### code chunk number 16: figE5
###################################################
plot(out, lwd = 2, ask = FALSE)
mtext(outer = TRUE, side = 3, line = -1.5, "E5", cex = 1.5)


###################################################
### code chunk number 17: deTestSet.Rnw:238-239
###################################################
out <- emep()


###################################################
### code chunk number 18: emep
###################################################
plot(out, which = 1:9, lwd = 2, ask = FALSE)
mtext(outer = TRUE, side = 3, line = -1.5, "emep", cex = 1.5)


###################################################
### code chunk number 19: figemep
###################################################
plot(out, which = 1:9, lwd = 2, ask = FALSE)
mtext(outer = TRUE, side = 3, line = -1.5, "emep", cex = 1.5)


###################################################
### code chunk number 20: deTestSet.Rnw:262-263
###################################################
out <- fekete()


###################################################
### code chunk number 21: fekete
###################################################
plot(out, which = 1:9, lwd = 2, ask = FALSE)
mtext(outer = TRUE, side = 3, line = -1.5, "fekete", cex = 1.5)


###################################################
### code chunk number 22: figfekete
###################################################
plot(out, which = 1:9, lwd = 2, ask = FALSE)
mtext(outer = TRUE, side = 3, line = -1.5, "fekete", cex = 1.5)


###################################################
### code chunk number 23: deTestSet.Rnw:286-287
###################################################
out <- hires()


###################################################
### code chunk number 24: hires
###################################################
plot(out, lwd = 2, ask = FALSE)
mtext(outer = TRUE, side = 3, line = -1.5, "hires", cex = 1.5)


###################################################
### code chunk number 25: fighires
###################################################
plot(out, lwd = 2, ask = FALSE)
mtext(outer = TRUE, side = 3, line = -1.5, "hires", cex = 1.5)


###################################################
### code chunk number 26: deTestSet.Rnw:310-311
###################################################
out <- nand(method = daspk)


###################################################
### code chunk number 27: nand
###################################################
plot(out, which = 1:9, lwd = 2, ask = FALSE)
mtext(outer = TRUE, side = 3, line = -1.5, "nand", cex = 1.5)


###################################################
### code chunk number 28: fignand
###################################################
plot(out, which = 1:9, lwd = 2, ask = FALSE)
mtext(outer = TRUE, side = 3, line = -1.5, "nand", cex = 1.5)


###################################################
### code chunk number 29: deTestSet.Rnw:334-335
###################################################
out <- orego()


###################################################
### code chunk number 30: orego
###################################################
plot(out, lwd = 2, ask = FALSE, log = "y")
mtext(outer = TRUE, side = 3, line = -1.5, "orego", cex = 1.5)


###################################################
### code chunk number 31: figorego
###################################################
plot(out, lwd = 2, ask = FALSE, log = "y")
mtext(outer = TRUE, side = 3, line = -1.5, "orego", cex = 1.5)


###################################################
### code chunk number 32: deTestSet.Rnw:358-359
###################################################
out <- pollution()


###################################################
### code chunk number 33: pollution
###################################################
plot(out, which = 1:9, lwd = 2, ask = FALSE)
mtext(outer = TRUE, side = 3, line = -1.5, "pollution", cex = 1.5)


###################################################
### code chunk number 34: figpollution
###################################################
plot(out, which = 1:9, lwd = 2, ask = FALSE)
mtext(outer = TRUE, side = 3, line = -1.5, "pollution", cex = 1.5)


###################################################
### code chunk number 35: deTestSet.Rnw:382-383
###################################################
out <- ring()


###################################################
### code chunk number 36: ring
###################################################
plot(out, which = 1:9, lwd = 2, ask = FALSE)
mtext(outer = TRUE, side = 3, line = -1.5, "ring", cex = 1.5)


###################################################
### code chunk number 37: figring
###################################################
plot(out, which = 1:9, lwd = 2, ask = FALSE)
mtext(outer = TRUE, side = 3, line = -1.5, "ring", cex = 1.5)


###################################################
### code chunk number 38: deTestSet.Rnw:406-407
###################################################
out <- rober()


###################################################
### code chunk number 39: rober
###################################################
plot(out, lwd = 2, ask = FALSE, log = "x", xlim = c(1e-5,1e11))
mtext(outer = TRUE, side = 3, line = -1.5, "rober", cex = 1.5)


###################################################
### code chunk number 40: figrober
###################################################
plot(out, lwd = 2, ask = FALSE, log = "x", xlim = c(1e-5,1e11))
mtext(outer = TRUE, side = 3, line = -1.5, "rober", cex = 1.5)


###################################################
### code chunk number 41: deTestSet.Rnw:430-431
###################################################
out <- transistor()


###################################################
### code chunk number 42: transistor
###################################################
plot(out, lwd = 2, ask = FALSE)
mtext(outer = TRUE, side = 3, line = -1.5, "transistor", cex = 1.5)


###################################################
### code chunk number 43: figtransistor
###################################################
plot(out, lwd = 2, ask = FALSE)
mtext(outer = TRUE, side = 3, line = -1.5, "transistor", cex = 1.5)


###################################################
### code chunk number 44: deTestSet.Rnw:454-455
###################################################
out <- tube()


###################################################
### code chunk number 45: tube
###################################################
plot(out, which = 1:9, lwd = 2, ask = FALSE)
mtext(outer = TRUE, side = 3, line = -1.5, "tube", cex = 1.5)


###################################################
### code chunk number 46: figtube
###################################################
plot(out, which = 1:9, lwd = 2, ask = FALSE)
mtext(outer = TRUE, side = 3, line = -1.5, "tube", cex = 1.5)


###################################################
### code chunk number 47: deTestSet.Rnw:478-479
###################################################
out <- twobit()


###################################################
### code chunk number 48: twobit
###################################################
plot(out, which = 1:9, lwd = 2, ask = FALSE)
mtext(outer = TRUE, side = 3, line = -1.5, "twobit", cex = 1.5)


###################################################
### code chunk number 49: figtwobit
###################################################
plot(out, which = 1:9, lwd = 2, ask = FALSE)
mtext(outer = TRUE, side = 3, line = -1.5, "twobit", cex = 1.5)


###################################################
### code chunk number 50: deTestSet.Rnw:502-503
###################################################
out <- vdpol()


###################################################
### code chunk number 51: vdpol
###################################################
plot(out, lwd = 2, ask = FALSE)
mtext(outer = TRUE, side = 3, line = -1.5, "vdpol", cex = 1.5)


###################################################
### code chunk number 52: figvdpol
###################################################
plot(out, lwd = 2, ask = FALSE)
mtext(outer = TRUE, side = 3, line = -1.5, "vdpol", cex = 1.5)


###################################################
### code chunk number 53: deTestSet.Rnw:527-528
###################################################
out <- wheelset()


###################################################
### code chunk number 54: wheelset
###################################################
plot(out, which = 1:6, lwd = 2, ask = FALSE)
mtext(outer = TRUE, side = 3, line = -1.5, "wheelset", cex = 1.5)


###################################################
### code chunk number 55: figwheel
###################################################
plot(out, which = 1:6, lwd = 2, ask = FALSE)
mtext(outer = TRUE, side = 3, line = -1.5, "wheelset", cex = 1.5)


