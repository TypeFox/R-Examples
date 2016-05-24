### R code from vignette source 'HAC.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
options(prompt = "R> ", continue = "+  ", width = 70, useFancyQuotes = FALSE)
library("HAC")


###################################################
### code chunk number 2: fourfully
###################################################
Obj1 = hac.full(type = 1, y = c("u4", "u3", "u2", "u1"), theta = c(2, 3, 4))
par(mai = c(0, 0, 0, 0))
plot(Obj1, index = TRUE, l = 1.6)


###################################################
### code chunk number 3: fourpartially
###################################################
Obj2 = hac(type = 1, tree = list(list("u4", "u3", 3), list("u1", "u2", 4), 2))
par(mai = c(0, 0, 0, 0))
plot(Obj2, index = TRUE, l = 1.6)


###################################################
### code chunk number 4: HAC.Rnw:402-406
###################################################
library("HAC")
data("finData")
system.time(result <- estimate.copula(finData, margins = "edf"))
result


###################################################
### code chunk number 5: Scatter1
###################################################
par(mai = c(0, 0, 0, 0))
pairs(finData, pch = 20, cex.axis = 1.5)


###################################################
### code chunk number 6: HAC.Rnw:515-516
###################################################
names(formals(estimate.copula))


###################################################
### code chunk number 7: result-agg
###################################################
    result.agg = estimate.copula(finData, method = 1, margins = "edf", epsilon = 0.3)
    par(mai = c(0, 0, 0, 0))
    plot(result.agg, circles = 0.3, index = TRUE, l = 1.7)


###################################################
### code chunk number 8: result
###################################################
    par(mai = c(0, 0, 0, 0))
    plot(result, circles = 0.3, index = TRUE, l = 1.7)


###################################################
### code chunk number 9: HAC.Rnw:602-605 (eval = FALSE)
###################################################
## result.agg = estimate.copula(sample, margins = "edf", epsilon = 0.3)
## plot(result, circles = 0.3, index = TRUE, l = 1.7)
## plot(result.agg, circles = 0.3, index = TRUE, l = 1.7)


###################################################
### code chunk number 10: HAC.Rnw:627-631
###################################################
G.cop = hac.full(type  = 1,
                 y     = c("X4", "X3", "X2", "X1"),
                 theta = c(1.1, 1.8, 2.5))
G.cop


###################################################
### code chunk number 11: HAC.Rnw:654-655
###################################################
hac(type = 1, tree = list("X1", "X2", "X3", "X4", 2))


###################################################
### code chunk number 12: HAC.Rnw:665-666
###################################################
hac(type = 1, tree = list(list("X1", "X2", 2.5), "X3", "X4", 1.5))


###################################################
### code chunk number 13: HAC.Rnw:676-680
###################################################
HAC = hac(type = 1, tree = list(list("Y1", list("Z3", "Z4", 3), "Y2", 2.5),
                      list("Z1", "Z2", 2), list("X1", "X2", 2.4),
                      "X3", "X4", 1.5))
HAC


###################################################
### code chunk number 14: HAC
###################################################
par(mai = c(0, 0, 0, 0))
plot(HAC, cex = 0.8, circles = 0.35)


###################################################
### code chunk number 15: HAC.Rnw:701-702
###################################################
plot(HAC, cex = 0.8, circles = 0.35)


###################################################
### code chunk number 16: HAC.Rnw:716-717
###################################################
names(formals(plot.hac))


###################################################
### code chunk number 17: Scatter2
###################################################
set.seed(1)
sim.data = rHAC(500, G.cop)
par(mai = c(0, 0, 0, 0))
pairs(sim.data, pch = 20, cex.axis = 1.75)


###################################################
### code chunk number 18: HAC.Rnw:791-793 (eval = FALSE)
###################################################
## sim.data = rHAC(500, G.cop)
## pairs(sim.data, pch = 20)


###################################################
### code chunk number 19: HAC.Rnw:818-819
###################################################
probs = pHAC(X = sim.data, hac = G.cop)


###################################################
### code chunk number 20: pp
###################################################
probs.emp = emp.copula.self(sim.data, proc = "M")
plot(probs, probs.emp, pch = 20, xlab = "True Probabilities", ylab = "Empirical Probabilites", asp = 1, cex.lab = 1.125, cex.axis = 1.125)
grid(lwd = 1.5)
box(lwd = 1)
points(probs, probs.emp, pch = 20)
lines(c(0,1), c(0,1), col = "darkred", lwd = 2)


###################################################
### code chunk number 21: HAC.Rnw:872-873 (eval = FALSE)
###################################################
## probs.emp = emp.copula.self(sim.data, proc = "M")


###################################################
### code chunk number 22: HAC.Rnw:879-883 (eval = FALSE)
###################################################
## emp.copula(u, x, proc = "M", sort = "none", margins = NULL,
##            na.rm = FALSE, ...)
## emp.copula.self(x, proc = "M", sort = "none", margins = NULL,
##                 na.rm = FALSE, ...)


