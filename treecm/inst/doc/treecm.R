### R code from vignette source 'treecm.Rnw'

###################################################
### code chunk number 1: ex1
###################################################
library(treecm)
data(stonePine1TreeData)
print(stonePine1TreeData)


###################################################
### code chunk number 2: ex11
###################################################
vectors  <- treeVectors(stonePine1TreeData)
CM       <- centreOfMass(vectors)
summary(CM)


###################################################
### code chunk number 3: ex2
###################################################
plot(vectors, main = "A stone pine centre of mass")
plot(CM)


###################################################
### code chunk number 4: ex3
###################################################
plot(vectors, 
  main = "A stone pine centre of mass",
  xlim = c(-8, 10),
  ylim = c(-12, 4)
)
plot(CM)
plotPolarSegment(210, 10.6, 140, 14.4)


###################################################
### code chunk number 5: exSnow
###################################################
rows <- substr(row.names(stonePine1TreeData$fieldData), 1, 1)
component <- substr(row.names(stonePine1TreeData$fieldData), 1, 1)
stonePine1TreeData$fieldData <- within(stonePine1TreeData$fieldData,
                                       biomass[height > 12 & component != "L"] <-  biomass[height > 12 & component != "L"] * 2
                                       )
rm("rows")


###################################################
### code chunk number 6: exSnow2
###################################################
vectors  <- treeVectors(stonePine1TreeData)
CM       <- centreOfMass(vectors)
summary(CM)
plot(vectors, 
  main = "A stone pine centre of mass under 2x snow load",
  xlim = c(-8,10),
  ylim = c(-12,4)
)
plot(CM)
plotPolarSegment(210, 10.6, 140, 14.4)


###################################################
### code chunk number 7: exWind
###################################################
data(stonePine1TreeData)

rows <- substr(row.names(stonePine1TreeData$fieldData), 1, 1)
stonePine1TreeData$fieldData <- within(
  stonePine1TreeData$fieldData, {
  biomass[((azimuth >= 270 | azimuth < 90) & rows != "L")] <- biomass[((azimuth >= 270 | azimuth < 90) & rows != "L")] / 2
  biomass[((azimuth >= 90 | azimuth < 270) & rows != "L")] <- biomass[((azimuth >= 90 | azimuth < 270) & rows != "L")] * 2  
})
rm(rows)

vectors  <- treeVectors(stonePine1TreeData)
CM       <- centreOfMass(vectors)
summary(CM)


###################################################
### code chunk number 8: exPruning
###################################################
library(treecm)
data(stonePine1TreeData)
vectors  <- treeVectors(stonePine1TreeData)
CM       <- centreOfMass(vectors)

op <- par(mfrow = c(2, 1), mai = c(0.5,0.5,0.5,0.2))

plot(vectors, main = "Centre of mass of a stone pine")
plot(CM)

component <- row.names(stonePine1TreeData$fieldData)
stonePine1TreeData$fieldData$toBePruned[component %in% c("B2", "B4")] <- TRUE
vectors  <- treeVectors(stonePine1TreeData)
CM       <- centreOfMass(vectors)
plot(vectors, main = "Centre of mass of a stone pine, B2 and B4 removed")
plot(CM)

par(op)
rm(op)


###################################################
### code chunk number 9: exSR
###################################################
data(stonePine1TreeData)
# assign length to branches
stonePine1TreeData$fieldData <- within(stonePine1TreeData$fieldData,
length[3:25] <- c(7, 7, 7, 7, 7, 7, 4, 7, 7, 4, 7, 7, 4, 7, 7, 7, 4, 7, 7, 4, 4, 7, 7)
)
vectors <- treeVectors(stonePine1TreeData)
SR      <- treeSR(stonePine1TreeData,vectors)
plot(SR, main = "Branches slenderness ratio", xaxt='n', yaxt = 'n', xlab = "", ylab = "")


###################################################
### code chunk number 10: exStabilization1
###################################################
library(treecm)
data(stonePine1TreeData)
vectors <- treeVectors(stonePine1TreeData)
CM      <- centreOfMass(vectors)

## We need to compute the tree moment
treeMoment <- buildTreeMomentObject(
  centreOfMassModulus(CM)
  , treeTotalBiomass(stonePine1TreeData)
  , centreOfMassAngle(CM)
)
treeMoment <- calcMoment(treeMoment)

## We extract the logs belonging to the main stem
mainStem <- logPathSelection(stonePine1TreeData)

(plinth <- getPlinthForce(
	l.stem = 10, 
	d = 40, 
	logs = mainStem, 
	treeMoment = getMoment(treeMoment), 
	CM = CM
))


###################################################
### code chunk number 11: exStabilization2
###################################################
plinth <- getPlinthForce(
	l.stem = 10, 
	d = 15:50, 
	logs = mainStem, 
	treeMoment = getMoment(treeMoment), 
	CM = CM
)
print(plinth$force)


###################################################
### code chunk number 12: exStabilization3
###################################################
plinth <- getPlinthForce(
	l.stem = seq(9, 12, 0.5), 
	d = 40, 
	logs = mainStem, 
	treeMoment = getMoment(treeMoment), 
	CM = CM
)
print(plinth$force)


###################################################
### code chunk number 13: exStabilization4
###################################################
aR <- anchorRange(mainStem, CM)
l.stemSeq <- round(seq(aR[["z"]] + 1, aR[["hMax"]] - 2, length.out = 6), 2)
plinth <- data.frame(
	getPlinthForce(
		l.stemSeq, 
		17:50, 
		mainStem, 
		getMoment(treeMoment), 
		CM
))
head(plinth)


###################################################
### code chunk number 14: exStabilization5
###################################################
library(ggplot2)
ggplot(data = plinth, aes(x = cableLength, y = force)) +
  geom_line(aes(color = factor(anchorAlongStem), group = anchorAlongStem)) + 
  ylab('Force [N]') +
  xlab('Cable length [m]') +
  labs(colour = "Anchor\nstem\ndistance [m]") +
  ggtitle("Force on the plinth (a stone pine)") +
	theme(
    legend.position = c(.8, .7), 
	  legend.background = element_rect(fill="white")
  )


###################################################
### code chunk number 15: exStabilization6
###################################################
plinth <- transform(plinth, 
                    distanceOverAnchorHeight = distanceOnGround / anchorAlongStem
                    , heightOverAnchorHeight = round(anchorAlongStem / CM[["z"]], 2)
                    , forceOverTreeBiomass   = force / treeTotalBiomass(stonePine1TreeData)/10
                    )
head(plinth[c("distanceOverAnchorHeight", "heightOverAnchorHeight", "forceOverTreeBiomass")])
ggplot(data = plinth, 
	aes(x = distanceOverAnchorHeight, y = forceOverTreeBiomass)) +
  geom_line() + 
  facet_wrap(~heightOverAnchorHeight)


###################################################
### code chunk number 16: dst
###################################################
library(treecm)
data(Dst)
print(Dst)


###################################################
### code chunk number 17: dstex1
###################################################
data(Dst)
with(Dst, density[species == "Pinus pinaster"])


###################################################
### code chunk number 18: allometryEx1
###################################################
csvFile <- system.file("doc", "CopyOfstonePine2FieldData.csv", package = "treecm")
sP.Cutini <- treeBiomass(importFieldData(csvFile, 650, allometryCutini2009))
head(sP.Cutini$fieldData)
sP.ABDC   <- treeBiomass(importFieldData(csvFile, 650, allometryABDC))
head(sP.ABDC$fieldData)

biomassRaw <- data.frame(
 cutini = sP.Cutini$fieldData$biomass
 , ABDC = sP.ABDC$fieldData$biomass
 , code = rownames(sP.ABDC$fieldData)
 , diameter = sP.ABDC$fieldData$dBase
 )
 
library(reshape2)
 
biomass <- melt(
 biomassRaw
 , measure.vars = c("ABDC", "cutini")
 , value.name = "biomass"
 , variable.name = "allometry"
)
 
rm(biomassRaw)
 
(treeBiomass <- with(biomass, tapply(biomass, allometry, sum)))


###################################################
### code chunk number 19: allometryBarPlot
###################################################
library(ggplot2)
ggplot(biomass, aes(x = code, y = biomass)) + 
  geom_bar(aes(fill = allometry), position = "dodge", stat = "identity")


###################################################
### code chunk number 20: allometryDotPlot
###################################################
library(ggplot2)
ggplot(biomass, aes(x = diameter, y = biomass)) + 
  geom_point(aes(shape = allometry), size = 5)


###################################################
### code chunk number 21: allometryEx2
###################################################
biomass <- within(biomass, {
  biomass[allometry == "cutini" & code != "L1"] <- biomass[allometry == "cutini" & code != "L1"] * 1.5
  biomass[allometry == "cutini" & code == "B4"] <- biomass[allometry == "ABDC" & code == "B4"]
})

with(biomass, tapply(biomass, allometry, sum))


###################################################
### code chunk number 22: allometryDotPlot2
###################################################
ggplot(biomass, aes(x = diameter, y = biomass)) + 
  geom_point(aes(shape = allometry), size = 5)


###################################################
### code chunk number 23: allometryEx3
###################################################
sP.Cutini$fieldData$biomass <- biomass$biomass[biomass$allometry == "cutini"]
head(sP.Cutini$fieldData)


