### R code from vignette source 'spTest-vig.rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: libs
###################################################
#install.packages("spTest")
#install.packages("geoR")
#install.packages("splines")
library("spTest")
library("fields")
library("geoR")
library("splines")


###################################################
### code chunk number 2: narccap
###################################################
data("WRFG")
image.plot(WRFG$lon-360, WRFG$lat, WRFG$WRFG.NCEP.tas, xlab = "Longitude",
ylab = "Latitude", main = "Mean WRFG-NCEP Temperatures")
world(add = T)


###################################################
### code chunk number 3: narccap2
###################################################
coords <- expand.grid(WRFG$xc,WRFG$yc)
sub <- which(coords[,1] > 2.95e6 & coords[,1] < 4.0e6 
& coords[,2] > 1.2e6 & coords[,2] < 2.25e6)
image.plot(WRFG$xc, WRFG$yc, WRFG$WRFG.NCEP.tas, xlab = "Easting", 
ylab = "Northing", main = "Map Projection Coordinates")
points(coords[sub,], pch = 16, cex = 0.35)
my.coords <- coords[sub,]
colnames(my.coords) <- c("xc", "yc") 
tas <- WRFG$WRFG.NCEP.tas
tas <- c(tas)
tas <- tas[sub]
mydata <- cbind(my.coords[,1],my.coords[,2], tas)


###################################################
### code chunk number 4: heatmap1
###################################################
x.coord <- unique(my.coords[,1])
y.coord <- unique(my.coords[,2])
nx <- length(x.coord)
ny <- length(y.coord)
tas.mat <- matrix(mydata[,3], nrow = nx, ncol = ny, byrow = F)
image.plot(x.coord, y.coord, tas.mat, ylab = "Northing", xlab = "Easting",
main = "Subset of Temperatures")


###################################################
### code chunk number 5: dirsemi1
###################################################
geodat <- as.geodata(mydata)
svar4 <- variog4(geodat)
plot(svar4)
title("Directional Sample Semivariograms")


###################################################
### code chunk number 6: test1
###################################################
my.delta <- 50000
mylags <-  rbind(c(4,0), c(0, 4), c(3, 3), c(-3,3))
myA <-  rbind(c(1, -1, 0 , 0), c(0, 0, 1, -1))
tr <- GuanTestGrid(spdata = mydata, delta = my.delta, lagmat = mylags, A = myA,
 window.dims = c(2,2),pt.est.edge = TRUE, sig.est.edge = TRUE, 
 sig.est.finite = TRUE, df = 2 )
tr$p.value.finite


###################################################
### code chunk number 7: model1
###################################################
x1 <- my.coords[,1]
x2 <- my.coords[,2]
m1 <- lm(tas ~ ns(x1, df = 3) + ns(x2, df = 3))
summary(m1)


###################################################
### code chunk number 8: iso2
###################################################
resid <- m1$resid
resid.mat <- matrix(resid, nrow = nx, ncol = ny, byrow = F)
image.plot(x.coord, y.coord, resid.mat, xlab = "Easting", ylab = "Northing")
title("Heat Map of Residuals")


###################################################
### code chunk number 9: dsvar2
###################################################
resid.dat <- cbind(mydata[,1:2], resid)
geodat.resid <- as.geodata(resid.dat)
plot(variog4(geodat.resid))
title("Directional Sample Semivariograms")


###################################################
### code chunk number 10: test2
###################################################
tr <- GuanTestGrid(spdata = resid.dat, delta = my.delta, lagmat = mylags, A = myA,
 window.dims = c(2,2),pt.est.edge = TRUE, sig.est.edge = TRUE, 
 sig.est.finite = TRUE, df = 2 )
tr$p.value.finite


###################################################
### code chunk number 11: rain
###################################################
data(COmonthlyMet)
sub30 <- CO.ppt[74:103,,]
nstations <- 376
years <- 1968:1997
nyears <- length(years)
yr.avg <- matrix(data = NA, nrow = nstations, ncol = nyears)
for(i in 1:nyears){
	yr.dat <- sub30[i,,]
	yr.avg[,i] <- apply(yr.dat, 2 , mean, na.rm = T)
}
avg30 <- apply(yr.avg, 1, mean, na.rm = T)
quilt.plot(CO.loc, log(avg30), xlab = "Longitude", ylab = "Latitude",
main = "Quilt Plot of log(precip)")
US(add = T)


###################################################
### code chunk number 12: elev
###################################################
quilt.plot(CO.loc, CO.elev, xlab = "Longitude", ylab = "Latitude",
main = "Quilt Plot of Elevation")
US(add = T)


###################################################
### code chunk number 13: rain1
###################################################
plot(CO.elev, log(avg30) , xlab = "Elevation", ylab = "log(precip)",
main = "Scatter of log(precip) vs. Elevation")
m1 <- lm( log(avg30) ~ ns(CO.elev, df = 3))
summary(m1)
fits <- m1$fitted.values
bad <- is.na(avg30)
x <- CO.elev[which(!bad)]
ord <- order(x)
x <- sort(x)
fits <- fits[ord]
lines(x, fits, lwd = 3, col = "red")


###################################################
### code chunk number 14: rain2
###################################################
resid <- m1$resid
hist(resid, freq = F, main = "Histogram of Residuals")
lines(density(resid), lwd = 2)


###################################################
### code chunk number 15: precipdv
###################################################
precip.resid <- cbind(CO.loc[which(!bad),][,1],CO.loc[which(!bad),][,2] , resid)
geodat <- as.geodata(precip.resid)
svar4 <- variog4(geodat)
plot(svar4, legend = F)
legend("bottomright", legend = c(expression(0*degree), expression(45*degree), expression(90*degree), expression(135*degree)), col = 1:4, lty = 1:4)
title("Directional Sample Semivariograms")


###################################################
### code chunk number 16: lags
###################################################
mylags <- 0.75 * rbind(c(1,0), c(0,1), c(1,1), c(-1,1))
myA <- rbind(c(1, -1, 0 , 0), c(0, 0, 1, -1))


###################################################
### code chunk number 17: grid
###################################################
quilt.plot(precip.resid[,1:2], precip.resid[,3], xlab = "Longitude", ylab = "Latitude")
title("Quilt Plot of Residuals and Grid Used for Subsampling")
min(precip.resid[,1])
max(precip.resid[,1])
min(precip.resid[,2])
max(precip.resid[,2])
#Set the limits of the sampling region to be slightly 
#larger than the min/max coordinates
my.xlims <-  c(-109.5, -101)
my.ylims <-  c(36.5, 41.5)
xlen <-  my.xlims[2]-my.xlims[1]
ylen <-  my.ylims[2]-my.ylims[1]
my.grid.spacing <-  c(xlen/12, ylen/12)
xgrid <- seq(my.xlims[1], my.xlims[2], by = my.grid.spacing[1])
ygrid <- seq(my.ylims[1], my.ylims[2], by = my.grid.spacing[2])
abline(v = xgrid, lty = 2)
abline(h = ygrid, lty = 2)


###################################################
### code chunk number 18: nongrid
###################################################
myh <-  0.6
myh.sb <-  0.7
tr.guan <- GuanTestUnif(spdata = precip.resid, lagmat = mylags, A = myA, df = 2, h = myh,
kernel = "norm", truncation = 1.5, xlims = my.xlims, ylims = my.ylims, 
grid.spacing = my.grid.spacing, window.dims = c(4,3), subblock.h = myh.sb)
tr.guan$p.value

tr.maity <- MaityTest(spdata = precip.resid, lagmat = mylags, A = myA, df = 2, 
xlims = my.xlims, ylims = my.ylims, grid.spacing = my.grid.spacing, 
block.dims = c(4,3), nBoot = 100, kappa = 1)
tr.maity$p.value


