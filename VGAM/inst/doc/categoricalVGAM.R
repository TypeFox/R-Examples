### R code from vignette source 'categoricalVGAM.Rnw'

###################################################
### code chunk number 1: categoricalVGAM.Rnw:84-90
###################################################
library("VGAM")
library("VGAMdata")
ps.options(pointsize = 12)
options(width = 72, digits = 4)
options(SweaveHooks = list(fig = function() par(las = 1)))
options(prompt = "R> ", continue = "+")


###################################################
### code chunk number 2: pneumocat
###################################################
pneumo <- transform(pneumo, let = log(exposure.time))
fit <- vgam(cbind(normal, mild, severe) ~ s(let, df = 2),
            cumulative(reverse = TRUE, parallel = TRUE), data = pneumo)


###################################################
### code chunk number 3: categoricalVGAM.Rnw:903-907
###################################################
journal <- c("Biometrika", "Comm.Statist", "JASA", "JRSS-B")
squaremat <- matrix(c(NA, 33, 320, 284,   730, NA, 813, 276,
                      498, 68, NA, 325,   221, 17, 142, NA), 4, 4)
dimnames(squaremat) <- list(winner = journal, loser = journal)


###################################################
### code chunk number 4: categoricalVGAM.Rnw:1007-1011
###################################################
abodat <- data.frame(A = 725, B = 258, AB = 72, O = 1073)
fit <- vglm(cbind(A, B, AB, O) ~ 1, ABO, data = abodat)
coef(fit, matrix = TRUE)
Coef(fit)  # Estimated pA and pB


###################################################
### code chunk number 5: categoricalVGAM.Rnw:1289-1291
###################################################
head(marital.nz, 4)
summary(marital.nz)


###################################################
### code chunk number 6: categoricalVGAM.Rnw:1294-1296
###################################################
fit.ms <- vgam(mstatus ~ s(age, df = 3), multinomial(refLevel = 2),
               data = marital.nz)


###################################################
### code chunk number 7: categoricalVGAM.Rnw:1300-1302
###################################################
head(depvar(fit.ms), 4)
colSums(depvar(fit.ms))


###################################################
### code chunk number 8: categoricalVGAM.Rnw:1311-1323
###################################################
# Plot output
mycol <- c("red", "darkgreen", "blue")
par(mfrow = c(2, 2))
plot(fit.ms, se = TRUE, scale = 12,
         lcol = mycol, scol = mycol)

# Plot output overlayed
#par(mfrow=c(1,1))
plot(fit.ms, se = TRUE, scale = 12,
         overlay = TRUE,
         llwd = 2,
         lcol = mycol, scol = mycol)


###################################################
### code chunk number 9: categoricalVGAM.Rnw:1366-1379
###################################################
getOption("SweaveHooks")[["fig"]]()
# Plot output
mycol <- c("red", "darkgreen", "blue")
 par(mfrow = c(2, 2))
 par(mar = c(4.2, 4.0, 1.2, 2.2) + 0.1)
plot(fit.ms, se = TRUE, scale = 12,
         lcol = mycol, scol = mycol)

# Plot output overlaid
#par(mfrow = c(1, 1))
plot(fit.ms, se = TRUE, scale = 12,
         overlay = TRUE,
         llwd = 2,
         lcol = mycol, scol = mycol)


###################################################
### code chunk number 10: categoricalVGAM.Rnw:1399-1400
###################################################
plot(fit.ms, deriv=1, lcol=mycol, scale=0.3)


###################################################
### code chunk number 11: categoricalVGAM.Rnw:1409-1413
###################################################
getOption("SweaveHooks")[["fig"]]()
# Plot output
 par(mfrow = c(1, 3))
 par(mar = c(4.5, 4.0, 0.2, 2.2) + 0.1)
plot(fit.ms, deriv = 1, lcol = mycol, scale = 0.3)


###################################################
### code chunk number 12: categoricalVGAM.Rnw:1436-1448
###################################################
foo <- function(x, elbow = 50)
  poly(pmin(x, elbow), 2)

clist <- list("(Intercept)" = diag(3),
             "poly(age, 2)" = rbind(1, 0, 0),
             "foo(age)"     = rbind(0, 1, 0),
             "age"          = rbind(0, 0, 1))
fit2.ms <-
    vglm(mstatus ~ poly(age, 2) + foo(age) + age,
         family = multinomial(refLevel = 2),
         constraints = clist,
         data = marital.nz)


###################################################
### code chunk number 13: categoricalVGAM.Rnw:1451-1452
###################################################
coef(fit2.ms, matrix = TRUE)


###################################################
### code chunk number 14: categoricalVGAM.Rnw:1456-1463
###################################################
par(mfrow = c(2, 2))
plotvgam(fit2.ms, se = TRUE, scale = 12,
         lcol = mycol[1], scol = mycol[1], which.term = 1)
plotvgam(fit2.ms, se = TRUE, scale = 12,
         lcol = mycol[2], scol=mycol[2], which.term = 2)
plotvgam(fit2.ms, se = TRUE, scale = 12,
         lcol = mycol[3], scol = mycol[3], which.term = 3)


###################################################
### code chunk number 15: categoricalVGAM.Rnw:1474-1483
###################################################
getOption("SweaveHooks")[["fig"]]()
# Plot output
par(mfrow=c(2,2))
 par(mar=c(4.5,4.0,1.2,2.2)+0.1)
plotvgam(fit2.ms, se = TRUE, scale = 12,
         lcol = mycol[1], scol = mycol[1], which.term = 1)
plotvgam(fit2.ms, se = TRUE, scale = 12,
         lcol = mycol[2], scol = mycol[2], which.term = 2)
plotvgam(fit2.ms, se = TRUE, scale = 12,
         lcol = mycol[3], scol = mycol[3], which.term = 3)


###################################################
### code chunk number 16: categoricalVGAM.Rnw:1501-1502
###################################################
deviance(fit.ms) - deviance(fit2.ms)


###################################################
### code chunk number 17: categoricalVGAM.Rnw:1508-1509
###################################################
(dfdiff <- df.residual(fit2.ms) - df.residual(fit.ms))


###################################################
### code chunk number 18: categoricalVGAM.Rnw:1512-1513
###################################################
pchisq(deviance(fit.ms) - deviance(fit2.ms), df = dfdiff, lower.tail = FALSE)


###################################################
### code chunk number 19: categoricalVGAM.Rnw:1526-1537
###################################################
ooo <- with(marital.nz, order(age))
with(marital.nz, matplot(age[ooo], fitted(fit.ms)[ooo, ],
     type = "l", las = 1, lwd = 2, ylim = 0:1,
     ylab = "Fitted probabilities",
     xlab = "Age",  # main="Marital status amongst NZ Male Europeans",
     col = c(mycol[1], "black", mycol[-1])))
legend(x = 52.5, y = 0.62,  # x="topright",
       col = c(mycol[1], "black", mycol[-1]),
       lty = 1:4,
       legend = colnames(fit.ms@y), lwd = 2)
abline(v = seq(10,90,by = 5), h = seq(0,1,by = 0.1), col = "gray", lty = "dashed")


###################################################
### code chunk number 20: categoricalVGAM.Rnw:1552-1565
###################################################
getOption("SweaveHooks")[["fig"]]()
 par(mfrow = c(1,1))
 par(mar = c(4.5,4.0,0.2,0.2)+0.1)
ooo <- with(marital.nz, order(age))
with(marital.nz, matplot(age[ooo], fitted(fit.ms)[ooo,],
     type = "l", las = 1, lwd = 2, ylim = 0:1,
     ylab = "Fitted probabilities",
     xlab = "Age",
     col = c(mycol[1], "black", mycol[-1])))
legend(x = 52.5, y = 0.62,
       col = c(mycol[1], "black", mycol[-1]),
       lty = 1:4,
       legend = colnames(fit.ms@y), lwd = 2.1)
abline(v = seq(10,90,by = 5), h = seq(0,1,by = 0.1), col = "gray", lty = "dashed")


###################################################
### code chunk number 21: categoricalVGAM.Rnw:1599-1603
###################################################
# Scale the variables? Yes; the Anderson (1984) paper did (see his Table 6).
head(backPain, 4)
summary(backPain)
backPain <- transform(backPain, sx1 = -scale(x1), sx2 = -scale(x2), sx3 = -scale(x3))


###################################################
### code chunk number 22: categoricalVGAM.Rnw:1607-1608
###################################################
bp.rrmlm1 <- rrvglm(pain ~ sx1 + sx2 + sx3, multinomial, data = backPain)


###################################################
### code chunk number 23: categoricalVGAM.Rnw:1611-1612
###################################################
Coef(bp.rrmlm1)


###################################################
### code chunk number 24: categoricalVGAM.Rnw:1640-1641
###################################################
set.seed(123)


###################################################
### code chunk number 25: categoricalVGAM.Rnw:1644-1646
###################################################
bp.rrmlm2 <- rrvglm(pain ~ sx1 + sx2 + sx3, multinomial, data = backPain, Rank = 2,
                   Corner = FALSE, Uncor = TRUE)


###################################################
### code chunk number 26: categoricalVGAM.Rnw:1654-1658
###################################################
biplot(bp.rrmlm2, Acol = "blue", Ccol = "darkgreen", scores = TRUE,
#      xlim = c(-1, 6), ylim = c(-1.2, 4),  # Use this if not scaled
       xlim = c(-4.5, 2.2), ylim = c(-2.2, 2.2),  # Use this if scaled
       chull = TRUE, clty = 2, ccol = "blue")


###################################################
### code chunk number 27: categoricalVGAM.Rnw:1690-1698
###################################################
getOption("SweaveHooks")[["fig"]]()
# Plot output
 par(mfrow=c(1,1))
 par(mar=c(4.5,4.0,0.2,2.2)+0.1)

biplot(bp.rrmlm2, Acol = "blue", Ccol = "darkgreen", scores = TRUE,
#      xlim = c(-1,6), ylim = c(-1.2,4),  # Use this if not scaled
       xlim = c(-4.5,2.2), ylim = c(-2.2, 2.2),  # Use this if scaled
       chull = TRUE, clty = 2, ccol = "blue")


###################################################
### code chunk number 28: categoricalVGAM.Rnw:1812-1813
###################################################
iam(NA, NA, M = 4, both = TRUE, diag = TRUE)


