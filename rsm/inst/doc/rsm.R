### R code from vignette source 'rsm.Rtex'

###################################################
### code chunk number 1: rsm.Rtex:102-104
###################################################
library("rsm")
ChemReact


###################################################
### code chunk number 2: rsm.Rtex:108-110
###################################################
CR1 <- coded.data(ChemReact1, x1 ~ (Time - 85)/5, x2 ~ (Temp - 175)/5)
CR1


###################################################
### code chunk number 3: rsm.Rtex:113-114
###################################################
as.data.frame(CR1)


###################################################
### code chunk number 4: rsm.Rtex:120-121 (eval = FALSE)
###################################################
## CR <- coded.data(ChemReact, x1 ~ (Time - 85)/5, x2 ~ (Temp - 175)/5)


###################################################
### code chunk number 5: rsm.Rtex:127-128
###################################################
code2val(data.frame(x1 = c(0.25, 0.5), x2 = c(-1.5, -0.5)), codings(CR1))


###################################################
### code chunk number 6: rsm.Rtex:135-137
###################################################
bbd(3, n0 = 2, coding =
  list(x1 ~ (Force - 20)/3, x2 ~ (Rate - 50)/10, x3 ~ Polish - 4))


###################################################
### code chunk number 7: rsm.Rtex:198-200
###################################################
ccd.pick(5, n.c = c(8, 16), blks.c = c(1, 2, 4),
  wbr.s = 1:2, restrict = "N<=65")


###################################################
### code chunk number 8: rsm.Rtex:212-214
###################################################
des1 <- ccd (y1 + y2 ~ A + B + C + D,
  generators = E ~ - A * B * C * D, n0 = c(6, 1))


###################################################
### code chunk number 9: rsm.Rtex:220-222
###################################################
des10 <- ccd( ~ A + B + C + D + E,
  blocks = Blk ~ c(A * B * C, C * D * E), n0 = c(2, 4))


###################################################
### code chunk number 10: varfcn
###################################################
par(mfrow=c(1,2))
varfcn(des10, ~ Blk + SO(A,B,C,D,E), dist = seq(0, 3, by=.1))
varfcn(des10, ~ Blk + SO(A,B,C,D,E), dist = seq(0, 3, by=.1), contour = TRUE)


###################################################
### code chunk number 11: rsm.Rtex:249-250
###################################################
ccd(2, n0 = c(1,1), inscribed=TRUE, randomize=FALSE)


###################################################
### code chunk number 12: rsm.Rtex:262-264
###################################################
CR1.rsm <- rsm(Yield ~ FO(x1, x2), data = CR1)
summary(CR1.rsm)


###################################################
### code chunk number 13: rsm.Rtex:270-272
###################################################
CR1.rsmi <- update(CR1.rsm, . ~ . + TWI(x1, x2))
summary(CR1.rsmi)


###################################################
### code chunk number 14: rsm.Rtex:278-279
###################################################
( CR2 <- djoin(CR1, ChemReact2) )


###################################################
### code chunk number 15: rsm.Rtex:284-286
###################################################
CR2.rsm <- rsm(Yield ~ Block + SO(x1, x2), data = CR2)
summary(CR2.rsm)


###################################################
### code chunk number 16: rsm.Rtex:293-295
###################################################
heli.rsm <- rsm(ave ~ block + SO(x1, x2, x3, x4), data = heli)
summary(heli.rsm)


###################################################
### code chunk number 17: fig6
###################################################
par(mfrow = c(2, 3))
contour(heli.rsm, ~ x1 + x2 + x3 + x4, image = TRUE,
  at = summary(heli.rsm)$canonical$xs)


###################################################
### code chunk number 18: rsm.Rtex:323-324
###################################################
steepest(CR1.rsm, dist = c(0, 0.5, 1))


###################################################
### code chunk number 19: rsm.Rtex:333-334
###################################################
canonical.path(heli.rsm, dist = seq(-5, 5, by = 0.5))


###################################################
### code chunk number 20: rsm.Rtex:340-343
###################################################
CO = as.coded.data(codata,  x1 ~ (Ethanol - 0.2)/0.1,  x2 ~ A.F.ratio - 15)
names(CO)[3] = "CO.conc"
head(CO)


###################################################
### code chunk number 21: rsm.Rtex:346-348
###################################################
CO.rsm = rsm(CO.conc ~ SO(x1,x2), data = CO)
canonical(CO.rsm)


###################################################
### code chunk number 22: rsm.Rtex:360-361
###################################################
canonical(CO.rsm, threshold = .2)


###################################################
### code chunk number 23: rr-code (eval = FALSE)
###################################################
## contour(CO.rsm, x2 ~ x1, bounds = list(x1=c(-16,2), x2=c(-2,16)), 
##         zlim=c(-100,100), col="gray", decode = FALSE)
## lines(c(-1,1,1,-1,-1), c(-1,-1,1,1,-1), col="green") # design region
## points(x2 ~ x1, data=canonical.path(CO.rsm), col="red", pch=1+6*(dist==0))
## points(x2 ~ x1, data=canonical.path(CO.rsm, threshold=.2), 
##         col="blue", pch=1+6*(dist==0))


###################################################
### code chunk number 24: rising-ridge
###################################################
par(mar=.1+c(4,4,0,0))
contour(CO.rsm, x2 ~ x1, bounds = list(x1=c(-16,2), x2=c(-2,16)), 
        zlim=c(-100,100), col="gray", decode = FALSE)
lines(c(-1,1,1,-1,-1), c(-1,-1,1,1,-1), col="green") # design region
points(x2 ~ x1, data=canonical.path(CO.rsm), col="red", pch=1+6*(dist==0))
points(x2 ~ x1, data=canonical.path(CO.rsm, threshold=.2), 
        col="blue", pch=1+6*(dist==0))


