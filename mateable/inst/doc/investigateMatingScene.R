## ----include=FALSE-------------------------------------------------------
knitr::opts_chunk$set(fig.width=7.2, fig.height=5)

## ------------------------------------------------------------------------
library(mateable)
packageDescription("mateable")$Version

## ------------------------------------------------------------------------
str(eelr2012)

## ------------------------------------------------------------------------
eelr <- makeScene(eelr2012, startCol = "firstDay", endCol = "lastDay",
                       xCol = "Ecoord", yCol = "Ncoord", idCol = "tagNo")

## ---- fig.show = 'hold'--------------------------------------------------
plotScene(eelr, "s")
plotScene(eelr, "s", pch = 1)

## ------------------------------------------------------------------------
ePair <- pairDist(eelr) # matrix of distances
hist(ePair, 40) # visualize histogram

## ------------------------------------------------------------------------
eKnn <- kNearNeighbors(eelr, 6) # 1,2,...,6 nearest potential mate
str(eKnn)
hist(eKnn[, 3], 40, main = "Distance to 3rd nearest neighbor (m)")

## ---- collapse = TRUE----------------------------------------------------
eSum <- matingSummary(eelr)
eSum[c("minX", "minY", "k1")] # by index
# by name
eSum$minX
eSum$minY
eSum$k1 # mean distance to 1st nearest potential mate

## ------------------------------------------------------------------------
eProx <- proximity(eelr, "maxPropSqrd")
eProx$pop
hist(eProx$ind$proximity, 30)


## ------------------------------------------------------------------------
plotPotential(eProx)

## ------------------------------------------------------------------------
focalPlants <- c(17217, 17202, 14582, 15114, 7614, 1509, 17002)
plotPotential(eProx, plotType = "net", sub.ids = focalPlants)

## ------------------------------------------------------------------------
plotScene(eelr, "t")
plotScene(eelr, "t", sub = 'all', drawQuartiles = FALSE, text.cex = 0.5)

## ------------------------------------------------------------------------

eOver <- overlap(eelr, compareToSelf = T) # matrix of days overlapping
hist(eOver, 40, main = "Histogram of days overlapping between pair") 

## ---- fig.show = 'hold'--------------------------------------------------
eRecep <- receptivityByDay(eelr) # T/F receptive on each day
str(eRecep) # matrix
dailyReceptivitySummary <- receptivityByDay(eelr, summary = TRUE)
dailyReceptivitySummary # a named integer vector
plot(as.Date(names(dailyReceptivitySummary)), dailyReceptivitySummary,
     xlab = 'date', ylab = 'count of receptive individuals')


## ------------------------------------------------------------------------
eSync <- synchrony(eelr, "overlap")
hist(eSync$ind[, 2], 30)
abline(v = eSync$pop, col ="red", lwd = 2)
abline(v = synchrony(eelr, "overlap", averageType = "median")$pop,
       col = "blue", lwd = 2)

## ------------------------------------------------------------------------
plotPotential(eSync, sub.ids = focalPlants)

## ------------------------------------------------------------------------
plotScene(eelr, c('s','t'), sub = focalPlants, N = 4, text.cex = 0.5)

## ------------------------------------------------------------------------
plot3DScene(eelr)
plot3DScene(eelr, pt.cex = .5, sub = focalPlants, N = 4)


## ---- collapse = TRUE, results='hold'------------------------------------
eSum <- matingSummary(eelr)
eSum[c("meanSD", "sdSD", "meanDur", "sdDur", "peak")] # index
eSum[c("minX", "minY", "k1")]
eSum$minX # by name
eSum$minY
eSum$k1

## ------------------------------------------------------------------------
# make scene based off eelr summary information
simScene <- simulateScene(size = nrow(eelr), meanSD = eSum$meanSD,
                          sdSD = eSum$sdSD, meanDur = eSum$meanDur,
                          sdDur = eSum$sdDur, xRange = c(eSum$minX, eSum$maxX),
                          yRange = c(eSum$minY, eSum$maxY))

## ------------------------------------------------------------------------
sProx <- proximity(simScene, "maxPropSqrd")
sSync <- synchrony(simScene, "augspurger")
sComp <- compatibility(simScene, "si_echinacea")

## ---- fig.height=4, fig.show = 'hold'------------------------------------
plotScene(simScene)
plot3DScene(simScene, sub = 'random')

plotPotential(sSync)
plotPotential(sProx)

## ------------------------------------------------------------------------
plotPotential(sComp, density = FALSE)
plotPotential(sComp, plotType = c('net','heat'), density = FALSE)

## ------------------------------------------------------------------------
plot3DPotential(list(sProx, sSync), subject = "ind")
plot3DPotential(list(sProx, sSync), subject = "ind", sample = 'all')
plot3DPotential(list(sProx, sSync), subject = "ind", sample = 'random')

plot3DPotential(list(sProx, sSync, sComp))
plot3DPotential(list(sProx, sSync, sComp), sub.ids = c(21,4))

## ---- collapse = TRUE, results='hold'------------------------------------

simulatedCoral <- simulateScene(10, sAlleles = 2)
plotScene(simulatedCoral)
plot3DScene(simulatedCoral)
plotPotential(compatibility(simulatedCoral, "dioecious"))

## ---- collapse = TRUE----------------------------------------------------
simulatedCoral
str(simulatedCoral)
summary(simulatedCoral)

## ------------------------------------------------------------------------
simScene1 <- simulateScene(size = nrow(eelr), meanSD = eSum$meanSD,
                          sdSD = eSum$sdSD, meanDur = eSum$meanDur,
                          sdDur = eSum$sdDur, xRange = c(eSum$minX, eSum$maxX),
                          yRange = c(eSum$minY, eSum$maxY))
simScene2<- simulateScene(size = 1.5*nrow(eelr), meanSD = eSum$meanSD + 365,
                          sdSD = eSum$sdSD, meanDur = eSum$meanDur,
                          sdDur = eSum$sdDur, xRange = c(eSum$minX, eSum$maxX),
                          yRange = c(eSum$minY, eSum$maxY))
simScene3 <- simulateScene(size = 0.8*nrow(eelr), meanSD = eSum$meanSD + 730,
                          sdSD = eSum$sdSD, meanDur = eSum$meanDur,
                          sdDur = eSum$sdDur, xRange = c(eSum$minX, eSum$maxX),
                          yRange = c(eSum$minY, eSum$maxY))

multiYearScene <- list('2012' = simScene1,'2013' = simScene2, '2014' = simScene3)

## ----fig.height = 8------------------------------------------------------
plotScene(multiYearScene,sub = c(1,6,12,18,13,24,55,45,60), text.cex = 0.8)

## ----fig.height = 8------------------------------------------------------
plot3DScene(multiYearScene, pt.cex = 1.2, sub = c(1,6,12,18,13,24,55,45,60))

## ------------------------------------------------------------------------
syncMulti <- synchrony(multiYearScene, method = 'augs')
proxMulti <- proximity(multiYearScene, method = 'maxPropSqrd')
compatMulti <- compatibility(multiYearScene, method = 'si_echinacea')

str(syncMulti) # a list of lists

## ---- fig.height = 8-----------------------------------------------------
plotPotential(syncMulti, sub.ids = c(1,6,12,18,13,24,55,44,60))

## ---- fig.height = 8-----------------------------------------------------
plotPotential(proxMulti,sub.ids = c(1,6,12,18,13,24,55,45,60))

## ---- fig.height = 8-----------------------------------------------------
plot3DPotential(list(syncMulti, proxMulti, compatMulti), subject = 'ind',
                pt.cex = 1, sub.ids = c(1,6,12,18,13,24,55,45,60))

