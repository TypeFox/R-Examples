### R code from vignette source 'OceanView.rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
options(prompt = " ")
options(continue = " ")
options(width=75)
library(OceanView)


###################################################
### code chunk number 2: OceanView.rnw:108-109
###################################################
head (WSnioz, n = 2)


###################################################
### code chunk number 3: OceanView.rnw:114-117
###################################################
NO3 <- db2cross(WSnioz, row = "SamplingDateTimeREAL", 
         col = "Station", val = "DataValue", 
         subset = (VariableName == "WNO3"), df.row = 5)


###################################################
### code chunk number 4: NO3
###################################################
image2D(NO3, resfac = 3)


###################################################
### code chunk number 5: NO3
###################################################
image2D(NO3, resfac = 3)


###################################################
### code chunk number 6: OceanView.rnw:141-143
###################################################
 head(WSnioz.table, n = 2)
 Msummary(WSnioz.table)


###################################################
### code chunk number 7: WSnioza
###################################################
 Mplot(WSnioz.table, subset = Station == 1, 
   select = c("WNO3", "WNO2"), xlab = "Daynr") 


###################################################
### code chunk number 8: WSnioza
###################################################
 Mplot(WSnioz.table, subset = Station == 1, 
   select = c("WNO3", "WNO2"), xlab = "Daynr") 


###################################################
### code chunk number 9: WSnioz
###################################################
 Mplot(Msplit(WSnioz.table, "Station", subset = Station %in% c(1, 13)) , 
   select = c("WNO3", "WNO2", "WNH4", "WO2"), lty = 1, lwd = 2,
   xlab = "Daynr", log = c("y", "y", "y", ""), 
   legend = list(x = "left", title = "Station")) 


###################################################
### code chunk number 10: WSnioz
###################################################
 Mplot(Msplit(WSnioz.table, "Station", subset = Station %in% c(1, 13)) , 
   select = c("WNO3", "WNO2", "WNH4", "WO2"), lty = 1, lwd = 2,
   xlab = "Daynr", log = c("y", "y", "y", ""), 
   legend = list(x = "left", title = "Station")) 


###################################################
### code chunk number 11: OceanView.rnw:200-201
###################################################
changeres(var = volcano, x = 1:nrow(volcano), y = 1:ncol(volcano), resfac = 0.1)


###################################################
### code chunk number 12: OceanView.rnw:204-206
###################################################
remap(var = volcano, x = 1:nrow(volcano), y = 1:ncol(volcano), 
  xto = c(1, 20, 40), yto = c(2, 5))


###################################################
### code chunk number 13: OceanView.rnw:209-211
###################################################
extract(volcano, x = 1:nrow(volcano), y = 1:ncol(volcano),
  xyto = cbind(c(2, 5), c(5, 10)))


###################################################
### code chunk number 14: flows
###################################################
 par(mfrow = c(2, 2))
 x  <- seq(-1, 1, by = 0.2)
 y  <- seq(-1, 1, by = 0.2)
 dx <- outer(x, y , function(x, y) -y)
 dy <- outer(x, y , function(x, y) x)

# velocity plot, different color for up/downward pointing arrows
 F <- quiver2D(u = dx, v = dy, x = x, y = y, colvar = dx > 0, 
     col = c("red", "blue"), colkey = FALSE, arr.max = 0.3, arr.min = 0.1)
 legend("topright", bg = "white", 
     legend = paste("max = ", format(F$speed.max, digits = 2))) 
 names(F)
 
 quiver2D(u = dx, v = dy, x = x, y = y, colvar = sqrt(dx^2 + dy^2), 
     arr.max = 0.1, arr.min = 0.1, clab = "speed")
     
# flow paths
 flowpath(u = dx, v = dy, x = x, y = y, numarr = 3, 
   startx = 0.1, starty = 0.1)
 flowpath(u = dx, v = dy, x = x, y = y, col = "red", numarr = 2, 
   startx = c(0.9, -0.9), starty = c(0.0, 0.0), add = TRUE)

# vectorplots
 u <- rnorm(10)
 v <- rnorm(10)
 x <- y <- 1 : 10
 vectorplot(u = u, v = v, x = x, y = y, clim = c(0, 3), 
   colvar = sqrt(u^2 + v^2), arr = TRUE)
 points(x, y)


###################################################
### code chunk number 15: flows
###################################################
 par(mfrow = c(2, 2))
 x  <- seq(-1, 1, by = 0.2)
 y  <- seq(-1, 1, by = 0.2)
 dx <- outer(x, y , function(x, y) -y)
 dy <- outer(x, y , function(x, y) x)

# velocity plot, different color for up/downward pointing arrows
 F <- quiver2D(u = dx, v = dy, x = x, y = y, colvar = dx > 0, 
     col = c("red", "blue"), colkey = FALSE, arr.max = 0.3, arr.min = 0.1)
 legend("topright", bg = "white", 
     legend = paste("max = ", format(F$speed.max, digits = 2))) 
 names(F)
 
 quiver2D(u = dx, v = dy, x = x, y = y, colvar = sqrt(dx^2 + dy^2), 
     arr.max = 0.1, arr.min = 0.1, clab = "speed")
     
# flow paths
 flowpath(u = dx, v = dy, x = x, y = y, numarr = 3, 
   startx = 0.1, starty = 0.1)
 flowpath(u = dx, v = dy, x = x, y = y, col = "red", numarr = 2, 
   startx = c(0.9, -0.9), starty = c(0.0, 0.0), add = TRUE)

# vectorplots
 u <- rnorm(10)
 v <- rnorm(10)
 x <- y <- 1 : 10
 vectorplot(u = u, v = v, x = x, y = y, clim = c(0, 3), 
   colvar = sqrt(u^2 + v^2), arr = TRUE)
 points(x, y)


###################################################
### code chunk number 16: OceanView.rnw:267-268
###################################################
dim(Ltrans)


###################################################
### code chunk number 17: Chesa
###################################################
image2D(Chesapeake$lon, Chesapeake$lat, z = Chesapeake$depth, 
  col = grey(seq(1, 0., length.out = 100)), main = "Ltrans",
  colkey = list(plot = FALSE))
  
scatter2D(x = Ltrans[,1,], y = Ltrans[,2,], colvar = Ltrans[,3,], 
  pch = ".", cex = 2, add = TRUE, clab = "depth, m")


###################################################
### code chunk number 18: Chesa
###################################################
image2D(Chesapeake$lon, Chesapeake$lat, z = Chesapeake$depth, 
  col = grey(seq(1, 0., length.out = 100)), main = "Ltrans",
  colkey = list(plot = FALSE))
  
scatter2D(x = Ltrans[,1,], y = Ltrans[,2,], colvar = Ltrans[,3,], 
  pch = ".", cex = 2, add = TRUE, clab = "depth, m")


###################################################
### code chunk number 19: Ltrans2D
###################################################
lon <- Chesapeake$lon
lat <- Chesapeake$lat
depth <- Chesapeake$depth

par(mfrow = c(2, 2)) 
for (i in seq(10, 106, length.out = 4)) 
   tracers2D(Ltrans[, 1, i], Ltrans[, 2, i],  
             colvar = Ltrans[ ,4, i], col = c("green", "orange"),
             pch = 16, cex = 0.5, 
             image = list(x = lon, y = lat, z = depth), colkey = FALSE,
             main = paste("time ", i))


###################################################
### code chunk number 20: Ltrans2D
###################################################
lon <- Chesapeake$lon
lat <- Chesapeake$lat
depth <- Chesapeake$depth

par(mfrow = c(2, 2)) 
for (i in seq(10, 106, length.out = 4)) 
   tracers2D(Ltrans[, 1, i], Ltrans[, 2, i],  
             colvar = Ltrans[ ,4, i], col = c("green", "orange"),
             pch = 16, cex = 0.5, 
             image = list(x = lon, y = lat, z = depth), colkey = FALSE,
             main = paste("time ", i))


###################################################
### code chunk number 21: Ltrans3D
###################################################
lon <- Chesapeake$lon
lat <- Chesapeake$lat
depth <- Chesapeake$depth

par(mfrow = c(1, 2), mar = c(0, 0, 2, 0)) 
for (i in c(20, 100)) 
tracers3D(Ltrans[, 1, i], Ltrans[, 2, i], Ltrans[, 3, i], 
          colvar = Ltrans[ ,4, i], col = c("green", "orange"),
          pch = 16, cex = 0.5, 
          surf = list(x = lon, y = lat, z = -depth, scale = FALSE, 
           expand = 0.02, colkey = FALSE, shade = 0.3, colvar = depth),
          colkey = FALSE, main = paste("time ", i))


###################################################
### code chunk number 22: Ltrans3D
###################################################
lon <- Chesapeake$lon
lat <- Chesapeake$lat
depth <- Chesapeake$depth

par(mfrow = c(1, 2), mar = c(0, 0, 2, 0)) 
for (i in c(20, 100)) 
tracers3D(Ltrans[, 1, i], Ltrans[, 2, i], Ltrans[, 3, i], 
          colvar = Ltrans[ ,4, i], col = c("green", "orange"),
          pch = 16, cex = 0.5, 
          surf = list(x = lon, y = lat, z = -depth, scale = FALSE, 
           expand = 0.02, colkey = FALSE, shade = 0.3, colvar = depth),
          colkey = FALSE, main = paste("time ", i))


