## ---- warning=F, message=FALSE-------------------------------------------
library(GeoLight)
data(hoopoe1)

head(hoopoe1)

## ------------------------------------------------------------------------
hoopoe1$datetime <- as.POSIXct(strptime(hoopoe1$datetime, format = "%Y-%m-%d %H:%M:%S", tz = "GMT"))
str(hoopoe1)

## ------------------------------------------------------------------------
twl <- twilightCalc(hoopoe1$datetime, hoopoe1$light, LightThreshold = 1.5, ask = F)
head(twl)

## ---- warning = F--------------------------------------------------------
data(calib2)
  calib2$tFirst <- as.POSIXct(calib2$tFirst, tz = "GMT")
  calib2$tSecond <- as.POSIXct(calib2$tSecond, tz = "GMT")


## Breeding location
lon.calib <- 8
lat.calib <- 47.01

## ---- fig.height=5, fig.width=7------------------------------------------
angle <- getElevation(calib2, known.coord = c(lon.calib, lat.calib), lnorm.pars = T)[1]

## ---- message=FALSE, fig.height=5, fig.width=7---------------------------
crds0 <- coord(twl, degElevation = angle, tol = 0)
crds1 <- coord(twl, degElevation = angle, tol = 0.13)

plot(twl[,1], crds0[,2], type = "o", pch = 16, col = "firebrick", 
     xlab = "Time", ylab = "Latitude")
points(twl[,1], crds1[,2], type = "o", pch = 16, col = "cornflowerblue")
abline(v = as.POSIXct("2008-09-21"), lty = 2)
legend("topleft", c("equinox", "tol = 0", "tol = 0.13"), pch = c(NA, 16, 16), lty = c(2,1,1), 
       col = c("black", "firebrick", "cornflowerblue"))

## ---- fig.height=10, fig.width=7-----------------------------------------
cL <- changeLight(twl$tFirst, twl$tSecond, type = twl$type, quantile=0.95, summary = F, days = 2)

## ---- fig.height=8, fig.width = 6----------------------------------------
siteMap(crds = crds1, site = cL$site, xlim = c(-12, 25), ylim = c(0, 50))

## ---- fig.height=10, fig.width=7-----------------------------------------
mS <- mergeSites(twl, site = cL$site, degElevation = angle, distThreshold = 300)

## ---- fig.width=6, fig.height = 10---------------------------------------
siteMap(crds1, mS$site, type = "cross", hull = F, lwd = 4, cex = 2,
        xlim = c(-12, 15), ylim = c(-30, 60))

arrows(mS$summary[-c(4,5),2], mS$summary[-c(4,5),5], mS$summary[-c(4,5),2], mS$summary[-c(4,5),7], 
       lty = 1, length = 0, lwd = 3)
arrows(mS$summary[-c(4,5),4], mS$summary[-c(4,5),3], mS$summary[-c(4,5),6], mS$summary[-c(4,5),3], 
       lty = 1, length = 0, lwd = 3)
points(mS$summary[-c(4,5),2:3], pch = 21, bg = "white", type = "b", lwd = 2, cex = 1.5)


arrows(mS$summary[c(4,5),2], mS$summary[c(4,5),5], mS$summary[c(4,5),2], mS$summary[c(4,5),7], 
       lty = 2, length = 0, lwd = 1)
arrows(mS$summary[c(4,5),2], mS$summary[c(4,5),5], mS$summary[c(4,5),2], mS$summary[c(4,5),7], 
       lty = 2, length = 0, lwd = 1)
points(mS$summary[c(4,5),2:3], pch = 21, bg = "grey90", type = "b", lwd = 1.5, cex = 1)

arrows(4.9, 46.01, -2.88, 40.46, length = 0.25, lwd = 2)


legend("bottomleft", pch = 21, 
       c("Positions based on 'mergeSites' optimization \nroutine (with 95% confidence intervals)."), 
       lty = 1)
legend("topleft", pch = 16, 
       c("Positions based on 'coord' function \nwith median and the 25 and 75 percent variation."), 
       lty = 1)

## ------------------------------------------------------------------------
schedule(twl$tFirst, twl$tSecond, site = mS$site)

