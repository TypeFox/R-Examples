### R code from vignette source 'ufc.rnw'

###################################################
### code chunk number 1: ufc.rnw:4-5
###################################################
source("../../scripts/functions.R")


###################################################
### code chunk number 2: ufc.rnw:27-28
###################################################
ufc.tree <- read.csv("../../data/ufc.csv")


###################################################
### code chunk number 3: ufc.rnw:33-34 (eval = FALSE)
###################################################
## str(ufc.tree)


###################################################
### code chunk number 4: ufc.rnw:36-37
###################################################
str(ufc.tree, vec.len = 1)


###################################################
### code chunk number 5: ufc.rnw:43-44
###################################################
show.cols.with.na(ufc.tree)


###################################################
### code chunk number 6: UFC-base
###################################################
names(ufc.tree) <- 
  c("point","tree","species","dbh.mm","ht.dm")
ufc.tree$dbh.cm <- ufc.tree$dbh.mm / 10
ufc.tree$ba.m2 <- ufc.tree$dbh.cm^2 / 40000 * pi
ufc.tree$height.m <- ufc.tree$ht.dm / 10


###################################################
### code chunk number 7: UFC-base
###################################################
height.hat <- ht.fvs.ni.m(ufc.tree$species, ufc.tree$dbh.cm)
missing.hts <- is.na(ufc.tree$height.m)
ufc.tree$height.m[missing.hts] <- height.hat[missing.hts]


###################################################
### code chunk number 8: UFC-base
###################################################
ufc.tree$vol.m3 <- 
    with(ufc.tree, vol.fvs.ni.m3(species, dbh.cm, height.m))


###################################################
### code chunk number 9: UFC-base
###################################################
ufc.baf.met <- 7
ufc.tree$tf.ha <- ufc.baf.met / ufc.tree$ba.m2
ufc.tree$vol.m3.ha <- ufc.tree$vol.m3 * ufc.tree$tf.ha

ufc.SyRS.data <- 
  aggregate(x = list(vol.m3.ha = ufc.tree$vol.m3.ha),
            by = list(point = ufc.tree$point),
            FUN = sum, 
            na.rm = TRUE)
ufc.SyRS.data$weight <- 1
str(ufc.SyRS.data)


###################################################
### code chunk number 10: UFC-base
###################################################
locations <- data.frame(point = 1:144, 
                        north.n = rep(c(12:1),12),
                        east.n = rep(c(1:12), rep(12,12)))

locations$north <- (locations$north.n - 0.5) * 134.11
locations$east <- (locations$east.n - 0.5) * 167.64


###################################################
### code chunk number 11: ufc.rnw:148-149
###################################################
ufc.SyRS.data <- merge(ufc.SyRS.data, locations)


###################################################
### code chunk number 12: fig-ufc
###################################################
opar <- par(las=1, pty="s")
plot(ufc.SyRS.data$east, ufc.SyRS.data$north, 
     type = "n", axes = F, 
     xlim = c(0,max(ufc.SyRS.data$east)+167.64/2),
     ylim = c(0,max(ufc.SyRS.data$north)+134.11/2), 
     xlab = "West-East (m)", ylab = "South-North (m)", 
     main = expression(paste("Units of 50", m^3, "/ha")))
axis(1); axis(2)
grayrange <- range(ufc.SyRS.data$vol.m3.ha)
text(formatC(ufc.SyRS.data$vol.m3.ha/50, 
             format = "f", digits = 0), 
     x = ufc.SyRS.data$east, 
     y = ufc.SyRS.data$north, cex=1.5,
     col = gray(1 - (ufc.SyRS.data$vol.m3.ha+200)/800))
par(opar)


###################################################
### code chunk number 13: ufc
###################################################
opar <- par(las=1, pty="s")
plot(ufc.SyRS.data$east, ufc.SyRS.data$north, 
     type = "n", axes = F, 
     xlim = c(0,max(ufc.SyRS.data$east)+167.64/2),
     ylim = c(0,max(ufc.SyRS.data$north)+134.11/2), 
     xlab = "West-East (m)", ylab = "South-North (m)", 
     main = expression(paste("Units of 50", m^3, "/ha")))
axis(1); axis(2)
grayrange <- range(ufc.SyRS.data$vol.m3.ha)
text(formatC(ufc.SyRS.data$vol.m3.ha/50, 
             format = "f", digits = 0), 
     x = ufc.SyRS.data$east, 
     y = ufc.SyRS.data$north, cex=1.5,
     col = gray(1 - (ufc.SyRS.data$vol.m3.ha+200)/800))
par(opar)


